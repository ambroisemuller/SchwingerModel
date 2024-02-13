#ifndef SCHWINGER_H
#define SCHWINGER_H

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include "lattice.h"

typedef enum
{
   LPFR,OMF2,OMF4,
   INTEGRATORS
} integrator_t;

/*data structure to store all the parameters of the algorithm*/
typedef struct {
   double eps;      /*integration step-size*/
   double tlength;  /*trajectory length*/
   int nstep;       /*leapfrog steps per trajectory*/
   int ntherm ;     /*number of thermalization steps*/
   int ntraj ;      /*number of trajectories after thermalization*/
   int subtraj ;    /*number of sub-trajectories per trajectory (i.e. trajectories between samples)*/
   int npf;         /*number of pseudofermion fields*/
   integrator_t integrator; /*MD integration scheme*/ 
   double lambda;   /*free parameter OMF2 integrator*/
   double mu[100];  /*number of twisted masses*/
   char measfile[128]; /*name of the measurement file*/
   double res_act,res_frc; /*residues of the solver for action and force*/
} hmc_params_t;

/*data structure to store all the parameters of the action*/
typedef struct {
   double beta;
   double kappa;
} act_params_t;

typedef struct {
   double complex s[2];
} spinor;

#endif
