/*******************************************************************************
 *
 * File gauge.c
 *
 * Copyright (C) 2017 Stefan Schaefer, 2021 Mattia Dalla Brida
 *
 * This software is distributed under the terms of the GNU General Public
 * License (GPL)
 *
 * HMC field update
 *
 * The externally accessible functions are
 *
 * double hmc(hmc_params_t *hmc_params, act_params_t *act_params)
 *   HMC update algorithm with algorithmic parameters passed in hmc_params
 *   and physics parameters in ac_params
 *   returns the acceptance rate
 *  
 *  Internal functions
 *   static double mom_action(double m[V][D])
 *     Returns the momentum contribution to the HMC Hamiltonian
 *     0.5*(m,m)
 *
 *   static double mom_heat(void)
 *     Initializes the momentum field and returns its action
 *
 *   static double hamiltonian(void)
 *      Computes the Hamiltonian of the HMC evolution, adding
 *      the various action terms
 *
 *   static void move_gauge(double eps)
 *      Implements the gauge field update step
 *         U -> e^{-eps*mom} U
 *
 *    static void move_mom(double eps)
 *      Computes the forces and updates the momenta 
 *         pi -> -eps F
 *    
 *    static omelyan4(double tau, int nstep)
 *       4th order integration of a trajectory of length tau with nstep steps
 *
 *    static double fermion_heat(void)
 *       Calls the heatbath for the various pseudofermion fields and accumulates
 *       and returns the action
 *
 *     static double start_hmc(void)
 *       Initializes all HMC fields and returns the value of the Hamiltonian
 *
 *******************************************************************************/

#include "update.h"
#include "linalg.h"
#include "random.h"
#include "measurements.h"
#include <time.h>

static double kappa,beta;
static double res_act,res_frc;

static int npf; 
static double *mus;
static double frc[V][D],mom[V][D];
static spinor **pf;


static double mom_action(double m[V][D])
{
   int i,mu;
   double act=0.0;

   for (i=0;i<V;i++)
      for (mu=0;mu<D;mu++)
         act+=m[i][mu]*m[i][mu];

   return act/2.0;
}


static double mom_heat(void)
{
   int i,mu;
   double norm=0.0;

   for(i=0;i<V;i++)
   {   
      for(mu=0;mu<D;mu++)
      {
         gauss_rand(1,mom[i]+mu);
         mom[i][mu]*=sqrt(2);
         norm+=mom[i][mu]*mom[i][mu];
      }
   }

   return 0.5*norm;
}


static void flip_mom(void)
{
   int i,mu;

   for(i=0;i<V;i++)
      for(mu=0;mu<D;mu++)
         mom[i][mu]=-mom[i][mu];
}


static double hamiltonian(void)
{
   int i;
   double act=0.0,sf,sg;

   act=mom_action(mom);

   if (npf>0)
   {   
      for (i=0;i<npf-1;i++)
      {
         sf=fermion_action2(pf[i],kappa,mus[i],mus[i+1],res_act);
         act+=sf;
      }

      sf=fermion_action1(pf[npf-1],kappa,mus[npf-1],res_act);
      act+=sf;
   }

   sg=gauge_action(beta);
   act+=sg;

   return act;
}


static void move_gauge(double eps)
{
   int i,mu;

   for (i=0;i<V;i++) 
      for (mu=0;mu<D;mu++) 
         gauge[i][mu]-=eps*mom[i][mu];
}

#ifdef MDINT_DBG

static void move_mom(double eps)
{
   int i,nu;
   double frc2[V][D];
   double nrm,nlk,cpu_time,cpu_time_tot;
   clock_t start,end;

   nlk=(double)(D)*(double)(V);
   cpu_time_tot=0.0;
   zero_frc(frc);

   if (npf>0)
   {
      for (i=0;i<npf-1;i++)
      {
         start=clock();
         fermion_force2(frc2,pf[i],kappa,mus[i],mus[i+1],1.0,res_frc);
         add_frc(frc,frc2);
         end=clock();

         cpu_time=((double)(end-start))/CLOCKS_PER_SEC;
         cpu_time_tot+=cpu_time;
         nrm=norm_frc(frc2);
         nrm=sqrt(nrm/nlk);
         printf("Force FRF_TM2, pf %d: ",i);
         printf("nrm = %.4e, eps = % .2e, nrm*|eps| = %.4e, "
                "time = %.2e sec\n",nrm,eps,nrm*fabs(eps),cpu_time);
      }

      start=clock();
      fermion_force1(frc2,pf[npf-1],kappa,mus[npf-1],1.0,res_frc);
      add_frc(frc,frc2);
      end=clock();

      cpu_time=((double)(end-start))/CLOCKS_PER_SEC;
      cpu_time_tot+=cpu_time;
      nrm=norm_frc(frc2);
      nrm=sqrt(nrm/nlk);
      printf("Force FRF_TM1, pf %d: ",npf-1);
      printf("nrm = %.4e, eps = % .2e, nrm*|eps| = %.4e, "
             "time = %.2e sec\n",nrm,eps,nrm*fabs(eps),cpu_time);

      nrm=norm_frc(frc);
      nrm=sqrt(nrm/nlk);
      printf("Force FRF_TOT:       ");
      printf("nrm = %.4e, eps = % .2e, nrm*|eps| = %.4e, "
             "time = %.2e sec\n",nrm,eps,nrm*fabs(eps),cpu_time_tot);
   }

   start=clock();
   gauge_force(frc2,beta,1.0);
   add_frc(frc,frc2);
   end=clock();

   cpu_time=((double)(end-start))/CLOCKS_PER_SEC;
   nrm=norm_frc(frc2);
   nrm=sqrt(nrm/nlk);
   printf("Force FRG:           ");
   printf("nrm = %.4e, eps = % .2e, nrm*|eps| = %.4e, "
          "time = %.2e sec\n",nrm,eps,nrm*fabs(eps),cpu_time);

   for (i=0;i<V;i++) 
      for (nu=0;nu<D;nu++) 
         mom[i][nu]-=eps*frc[i][nu];
}

#else

static void move_mom(double eps)
{
   int i,nu;
   double frc2[V][D];

   zero_frc(frc);

   if (npf>0)
   {   
      for (i=0;i<npf-1;i++)
      {
         fermion_force2(frc2,pf[i],kappa,mus[i],mus[i+1],1.0,res_frc);
         add_frc(frc,frc2);
      }

      fermion_force1(frc2,pf[npf-1],kappa,mus[npf-1],1.0,res_frc);
      add_frc(frc,frc2);
   }

   gauge_force(frc2,beta,1.0);
   add_frc(frc,frc2);

   for (i=0;i<V;i++) 
      for (nu=0;nu<D;nu++) 
         mom[i][nu]-=eps*frc[i][nu];
}

#endif

static void omelyan4(double tau,int nstep)
{
   int i;
   double eps=tau/nstep;
   double r1=0.08398315262876693,
          r2=0.2539785108410595,
          r3=0.6822365335719091,
          r4=-0.03230286765269967;

   for (i=0;i<nstep;i++)
   {
      move_mom(r1*eps);
      move_gauge(r2*eps);
      move_mom(r3*eps);
      move_gauge(r4*eps);
      move_mom((0.5-r1-r3)*eps);
      move_gauge((1-2*(r2+r4))*eps);
      move_mom((0.5-r1-r3)*eps);
      move_gauge(r4*eps);
      move_mom(r3*eps);
      move_gauge(r2*eps);
      move_mom(r1*eps);
   }
}


static void omelyan2(double lambda,double tau,int nstep)
{
   int i;
   double eps=tau/nstep;

   for (i=0;i<nstep;i++)
   {
      move_mom(lambda*eps);
      move_gauge(0.5*eps);
      move_mom((1.0-2.0*lambda)*eps);
      move_gauge(0.5*eps);
      move_mom(lambda*eps);
   }
}


static void leapfrog(double tau,int nstep)
{
   int i;
   double eps=tau/nstep;

   for (i=0;i<nstep;i++)
   {
      move_mom(0.5*eps);
      move_gauge(eps);
      move_mom(0.5*eps);
   }
}


static double fermion_heat(void)
{
   int ipf;
   double act=0.0,c;

   if (npf>0)
   {
      for (ipf=0;ipf<npf-1;ipf++)
      {
         c=setpf2(pf[ipf],kappa,mus[ipf],mus[ipf+1],res_act);
         act+=c;
      }

      c=setpf1(pf[npf-1],kappa,mus[npf-1]);
      act+=c;
   }   

   return act;
}


static double start_hmc()
{
   double H;
   
   H=mom_heat();
   H+=fermion_heat();
   H+=gauge_action(beta);

   return H;
}


void reversibility(hmc_params_t *hmc_params,act_params_t *act_params)
{
   int i,nstep;
   double tau,lambda;
   double H1,H2,dH,ddH,Dg,Dp;
   double gauge_old[V][D],mom_old[V][D];

   tau=hmc_params->tlength;
   nstep=hmc_params->nstep;
   lambda=hmc_params->lambda;

   kappa=act_params->kappa;
   beta=act_params->beta;
   npf=hmc_params->npf;
   mus=hmc_params->mu;
   res_act=hmc_params->res_act;
   res_frc=hmc_params->res_frc;

   if (npf>0)
   {
      pf=malloc(npf*sizeof(spinor*));
      for (i=0;i<npf;i++)
         pf[i]=malloc(V*sizeof(spinor));
   }
   else
      pf=NULL;

   H1=start_hmc();
   assign_link(gauge,gauge_old);
   assign_link(mom,mom_old);

   if (hmc_params->integrator==LPFR)
      leapfrog(tau,nstep);
   else if (hmc_params->integrator==OMF2)
      omelyan2(lambda,tau,nstep);
   else if (hmc_params->integrator==OMF4)
      omelyan4(tau,nstep);
   else
      printf("Unknown integrator\n");

   H2=hamiltonian();
   dH=H2-H1;

   flip_mom();

   if (hmc_params->integrator==LPFR)
      leapfrog(tau,nstep);
   else if (hmc_params->integrator==OMF2)
      omelyan2(lambda,tau,nstep);
   else if (hmc_params->integrator==OMF4)
      omelyan4(tau,nstep);
   else
      printf("Unknown integrator\n");

   H2=hamiltonian();
   ddH=fabs(H2-H1);

   flip_mom();

   Dg=cmp_link(gauge,gauge_old);
   Dp=cmp_link(mom,mom_old);

   printf("dH %+4.3e   ddH %4.3e   Dg %15.12e   Dp %15.12e\n",dH,ddH,Dg,Dp);
   fflush(stdout);

   if (npf>0)
   {
      for (i=0;i<npf;i++)
         free(pf[i]);
      free(pf);
   }
}   


double hmc(hmc_params_t *hmc_params,act_params_t *act_params)
{
   int itraj,isubtraj,accept,acc,i,nstep;
   double H1,H2,dH,r;
   double tau,lambda;
   double gauge_old[V][D];

   tau=hmc_params->tlength/(double)hmc_params->subtraj;
   nstep=hmc_params->nstep;
   lambda=hmc_params->lambda;

   kappa=act_params->kappa;
   beta=act_params->beta;
   npf=hmc_params->npf;
   mus=hmc_params->mu;
   res_act=hmc_params->res_act;
   res_frc=hmc_params->res_frc;

   if (npf>0)
   {
      pf=malloc(npf*sizeof(spinor*));
      for (i=0;i<npf;i++)
         pf[i]=malloc(V*sizeof(spinor));
   }
   else
      pf=NULL;

   accept=0;

   datafile_headers(hmc_params, act_params);

   printf("tau %f \n", tau);

   for (itraj=0;itraj<hmc_params->ntraj;itraj++)
   {

      for (isubtraj=0; isubtraj<hmc_params->subtraj;isubtraj++)
      {
         assign_link(gauge,gauge_old);

         /* initialize pseudofermions and momenta */
         H1=start_hmc();

         /* MD integration */
         if (hmc_params->integrator==LPFR)
            leapfrog(tau,nstep);
         else if (hmc_params->integrator==OMF2)
            omelyan2(lambda,tau,nstep);
         else if (hmc_params->integrator==OMF4)
            omelyan4(tau,nstep);
         else
            printf("Unknown integrator\n");

         /* acceptance step */
         H2=hamiltonian();
         dH=H2-H1;
         acc=0;

         if (dH<0)
            acc=1;
         else
         {
            ranlxd(&r,1);
            if (exp(-dH)>r)
               acc=1;
            else
               assign_link(gauge_old,gauge);
         }  
         printf("traj %d  subtraj %d  dH %+4.3e  acc %d\n", itraj, isubtraj, dH, acc);
         fflush(stdout);
         accept+=acc;
      }

      printf("measuring correlators...\n");
      measure_and_record(pf, itraj*tau*(double)hmc_params->subtraj, dH, acc);

      fflush(stdout);
      
   }

   if (npf>0)
   {
      for (i=0;i<npf;i++)
         free(pf[i]);
      free(pf);
   }

   run_plot_scripts();

   return accept/((double)hmc_params->ntraj*(double)hmc_params->subtraj);
}

