
/* 
 *   File schwinger.c
 *
 *   Contains the main program and a few other routines to simulate the 
 *   Schwinger model
 *   Routines for reading in the main parameters of the action and the 
 *   algorithm as  well as the computation of the action are provided.
 *
 *   This code is run with either one or two arguments
 *
 */ 

#define CONTROL
#include "lattice.h"
#include "measurements.h"
#include "schwinger.h"
#include "geometry.h"
#include "random.h"
#include "update.h"
#include <string.h>

/*  
 *  data structures to store all the parameters of the algorithm,
*/

static hmc_params_t hmc_params;
static act_params_t act_params;
static int seed;
char cfile[128];


static int get_val(FILE* fp, char *str, char* fmt,  void* val)
{
   char c[128];

   if(1!=fscanf(fp,"%s",c))
   {
      fprintf(stderr,"Error reading input file at %s\n",str);
      exit(1);
   }

   if(strcmp(str,c)!=0)
   {
      fprintf(stderr,"Error reading input file expected %s found %s\n",str,c);
      exit(1);
   }

   if(1!=fscanf(fp,fmt,val))
   {
      fprintf(stderr,"Error reading input file at %s\n",str);
      fprintf(stderr,"Cannot read value format %s\n",fmt);
      exit(1);
   }

   return 0;
}


static int read_input(char *input)
{
   int i;
   char c[128];
   FILE *fp;

   fp=fopen(input,"r");

   if (fp==NULL) 
   {
      fprintf(stderr, "Cannot open input file %s \n",input);
      exit(1);
   }

   get_val(fp, "kappa",       "%lf",&act_params.kappa);
   get_val(fp, "beta",        "%lf",&act_params.beta);
   get_val(fp, "ntherm",      "%i",&hmc_params.ntherm);
   get_val(fp, "ntraj",       "%i",&hmc_params.ntraj);
   get_val(fp, "subtraj",       "%i",&hmc_params.subtraj);
   get_val(fp, "traj_length", "%lf",&hmc_params.tlength);
   get_val(fp, "nstep",       "%i",&hmc_params.nstep);
   get_val(fp, "integrator",  "%s",c);
   if (strcmp(c,"LPFR")==0)
   {
      hmc_params.integrator=LPFR;
      hmc_params.lambda=0.0;
   }
   else if (strcmp(c,"OMF2")==0)
   {
      hmc_params.integrator=OMF2;
      get_val(fp,"lambda",    "%lf",&hmc_params.lambda);
   }
   else if (strcmp(c,"OMF4")==0)
   {
      hmc_params.integrator=OMF4;
      hmc_params.lambda=0.0;
   }
   else
   {
      fprintf(stderr,"Unknown integrator %s",c);
      exit(1);
   }
   get_val(fp, "npf", "%i",&hmc_params.npf);
   for (i=0;i<hmc_params.npf;i++)
      get_val(fp, "mu", "%lf",hmc_params.mu+i);
   get_val(fp, "res_act", "%lf",&hmc_params.res_act);
   get_val(fp, "res_frc", "%lf",&hmc_params.res_frc);
   get_val(fp, "seed",        "%i",&seed);
   get_val(fp, "config",    "%s",cfile);

   printf("PARAMETERS\n");
   printf("L              %i\n",L);
   printf("T              %i\n",T);
   printf("DIM            %i\n",D);
   printf("kappa          %f\n",act_params.kappa);
   printf("beta           %f\n",act_params.beta);
   printf("ntherm         %i\n",hmc_params.ntherm);
   printf("ntraj          %i\n",hmc_params.ntraj);
   printf("subtraj          %i\n",hmc_params.subtraj);
   printf("traj_length    %f\n",hmc_params.tlength);
   if (hmc_params.integrator==LPFR)
      printf("Leapfrog integrator\n");
   else if (hmc_params.integrator==OMF2)
      printf("2nd order OMF integrator with lambda = %e\n",
             hmc_params.lambda);
   else if (hmc_params.integrator==OMF4)
      printf("4th order OMF integrator\n");
   else
      printf("Unknown integrator\n");
   printf("nstep          %i\n",hmc_params.nstep);
   printf("npf            %i\n",hmc_params.npf);
   for (i=0;i<hmc_params.npf;i++)
      printf("mu             %f\n",hmc_params.mu[i]);
   printf("res_act        %e\n",hmc_params.res_act);
   printf("res_frc        %e\n",hmc_params.res_frc);
   printf("seed           %d\n",seed);
   printf("config_file    %s\n",cfile);
   printf("END PARAMETERS\n");
   fclose(fp);

   return 0;
}


static void save_config(char *cfile)
{
   int nw;
   FILE *fp;

   fp=fopen(cfile,"wb");

   if (fp==NULL) 
   {
      fprintf(stderr,"Error writing config");
      exit(1);
   }

   nw=fwrite(gauge,sizeof(double),V*D,fp);
   fclose(fp);

   if (nw!=V*D) 
   {
      fprintf(stderr,"Error writing config");
      exit(1);
   }
}


static void read_config(char *cfile)
{
   int nw;
   FILE *fp;

   fp=fopen(cfile,"rb");
   nw=fread(gauge,sizeof(double),V*D,fp);
   fclose(fp);

   if (nw!=V*D) 
   {
      fprintf(stderr,"Error writing config");
      exit(1);
   }
}


int main(int argc, char* argv[])
{
   double acc;
   int ix, mu;
   char config_file[128];
   char bc[2];

   if (PERIODIC == 0){
        bc[0] = 'A';
    } else {
        bc[0] = 'P';
    }
   bc[1] = '\0';

   if (argc!=2 && argc!=3) 
   {
      fprintf(stderr, "Number of arguments not correct\n");
      fprintf(stderr, "Usage: %s <infile> \n",argv[0]);
      exit(1);
   }

   /* get the parameters from the input file */
   read_input(argv[1]);

   /* initialize random number generator */
   rlxd_init(1,seed);

   /* initialize the nearest neighbor field */
   hopping(hop);

   /* initialize gauge field, either cold start or read in config */
   if (argc==2)
   {
      for (ix=0;ix<V;ix++)
      {   
         for (mu=0;mu<D;mu++)
            gauge[ix][mu]=0;
      } 

      printf("COLD START\n");
   }    
   else
   {
      sprintf(config_file, "%sT%d_L%d_%s/k%.4f_b%.3f/%s", DATA_FOLDER, T, L, bc, act_params.kappa, act_params.beta, argv[2]);
      read_config(config_file);
      printf("READ CONF. %s\n",config_file);
   }

   /* do hmc */
   acc=hmc(&hmc_params,&act_params);
   printf("Acc-rate %4.3e\n",acc);

   sprintf(config_file, "%sT%d_L%d_%s/k%.4f_b%.3f/%s", DATA_FOLDER, T, L, bc, act_params.kappa, act_params.beta, cfile);
   save_config(config_file);
   printf("\nField configurations saved to %s\n", config_file);

   return 0;
}

