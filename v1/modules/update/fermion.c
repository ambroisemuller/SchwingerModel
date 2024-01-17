/*******************************************************************************
 *
 * File fermion.c
 *
 * Copyright (C) 2017 Stefan Schaefer
 *
 * This software is distributed under the terms of the GNU General Public
 * License (GPL)
 *
 * For the 2D Schwinger model: Wilons fermion action and forces
 * with (twisted mass) Hasenbusch ratios
 *
 * The externally accessible functions are
 *
 * double setpf1(spinor *pf,double kappa,double mu)
 *   Fermion heatbath for pseudofermion field pf, corresponding to 
 *   fermion_action1()
 *
 * void fermion_force1(double f[V][D],spinor *phi,double kappa,double mu,double a)
 *   Sets the force field f to a times the derivative of the twisted mass 
 *   fermion_action1() for pseudofermion field phi with mass parameter 
 *   kappa and twisted mass mu
 *
 * double fermion_action1(spinor *phi,double kappa,double mu)
 *   Fermion action Sf1=(phi,Q^{-2} phi)
 *
 * double setpf2(spinor *pf,double kappa,double mu1,double mu2)
 *   Heatbath for the ratio of fermion matrices, corresponding to fermion_action2()
 *
 * void fermion_force2(double f[V][D],spinor *phi,double kappa,double mu1,
 *                     double mu2,double a)
 *   Sets the force field f to a times the derivative of the twisted mass 
 *   fermion_action2() for pseudofermion field phi with mass parameter 
 *   kappa and twisted mass mu
 *   
 * double fermion_action2(spinor *phi,double kappa,double mu1,double mu2)
 *   Action for determinant ratio Sf2=(phi,(Q^2+mu2^2)/(Q^2+mu1^2) phi)
 *   
 *******************************************************************************/

#include "update.h"
#include "linalg.h"
#include "dirac.h"
#include "random.h"

/* workspace */
static spinor v1[V],v2[V];


double setpf1(spinor *pf,double kappa,double mu)
{
   double norm;

   gauss_spinor(v1);
   dirac(v1,pf,kappa,mu);
   gamma5(pf);
   norm=creal(scalar_prod(v1,v1));

   return norm;
}


void fermion_force1(double f[V][D],spinor *phi,double kappa,double mu,double a,double res)
{
   int ix,nu;
   spinor t;
   complex double c;

   /* v1=Q^-2 phi */
   cg(phi,v1,kappa,mu,res,1000);

   /* v2=Q^-1 phi */
   dirac(v1,v2,kappa,mu);
   gamma5(v2);

   for (ix=0;ix<V;ix++)
   {
      for (nu=0;nu<D;nu++)
      {
         c=-a*kappa*cexp(-I*gauge[ix][nu]);
         mul_1_plus_gamma(nu,v2+hop[ix][nu],&t);
         f[ix][nu]=2*(cimag(conj(v1[ix].s[0])*c*t.s[0])-cimag(conj(v1[ix].s[1])*c*t.s[1]));

         mul_1_plus_gamma(nu,v1+hop[ix][nu],&t);
         f[ix][nu]+=2*(cimag(conj(v2[ix].s[0])*c*t.s[0])-cimag(conj(v2[ix].s[1])*c*t.s[1]));
      }
   }
}


double fermion_action1(spinor *phi,double kappa,double mu,double res)
{
   cg(phi,v1,kappa,mu,res,1000);

   return creal(scalar_prod(phi,v1));
}


double setpf2(spinor *pf,double kappa,double mu1,double mu2,double res)
{
   double norm;

   gauss_spinor(pf);

   cg(pf,v1,kappa,mu2,res,1000);
   assign_spinor(v1,v2);
   dirac(v2,v1,kappa,-mu2);
   gamma5(v1);

   mulc_add_assign(pf,v1,I*(mu1-mu2));
   norm=(mu2*mu2-mu1*mu1)*creal(scalar_prod(v1,v1));

   return norm;
}


void fermion_force2(double f[V][D],spinor *phi,double kappa,double mu1,double mu2,double a,double res)
{
   fermion_force1(f,phi,kappa,mu1,a*(mu2*mu2-mu1*mu1),res);
}


double fermion_action2(spinor *phi,double kappa,double mu1,double mu2,double res)
{
   return (mu2*mu2-mu1*mu1)*fermion_action1(phi,kappa,mu1,res);
}

