/*******************************************************************************
 *
 * File dirac.c
 *
 * Copyright (C) 2017 Stefan Schaefer
 *
 * This software is distributed under the terms of the GNU General Public
 * License (GPL)
 *
 * For the 2D Schwinger model: Wilons Dirac operator and some helper functions
 * 
 * The externally accessible functions are
 *
 * void gamma5(spinor *s)
 *   Replace volume spinor field s with gamma5*s
 *
 * void mul_1_plus_gamma(int mu,spinor *s, spinor *r)
 *   r = (1+gamma_mu) s
 *
 * void dirac(spinor *s,spinor *r,double kappa,double mu)
 *   Dirac operator 
 *   On exit s=D*r
 *   Mass parameter kappa, twisted mass mu
 *
 *******************************************************************************/

#include "dirac.h"

#define PI 3.14159


void gamma5(spinor *s)
{
   int i;

   for (i=0;i<V;i++)
      s[i].s[1]=-s[i].s[1];
}


void mul_1_plus_gamma(int mu,spinor *s,spinor *r)
{
   if (mu==0)
   {
      (*r).s[0]=((*s).s[0]-(*s).s[1]);
      (*r).s[1]=-(*r).s[0];
   }
   else if (mu==1)
   {
      (*r).s[0]=((*s).s[0]+I*(*s).s[1]);
      (*r).s[1]=-I*(*r).s[0];
   }
   else if (mu==2)
   {
      (*r).s[0]=((*s).s[0]+(*s).s[1]);
      (*r).s[1]=(*r).s[0];
   } 
   else if (mu==3)
   {
      (*r).s[0]=((*s).s[0]-I*(*s).s[1]);
      (*r).s[1]=I*(*r).s[0];
   }
}


void dirac(spinor *s,spinor *r,double kappa,double mu)
{
   int ix, nu;
   double complex c;
   spinor t;
   double phase = 0;

   for (ix=0;ix<V;ix++)
   {
      r[ix].s[0]=(1+I*mu)*s[ix].s[0];
      r[ix].s[1]=(1-I*mu)*s[ix].s[1];

      for (nu=0;nu<2*D;nu++)
      {
         if (nu<D) /* forward hopping*/
         {
            /*
            if (nu == 0 && (ix/L == T-1)){
               phase = PI;
            } else{
               phase = 0;
            }
            */
            c=-kappa*cexp(-I*(gauge[ix][nu] + phase)); /* pi phase for time boundaries (both forward and backward) */
         }
         else     /*backward hopping*/
         {
            /*
            if (nu == D && (ix/L == 0)){
               phase = PI;
            } else{
               phase = 0;
            }
            */
            c=-kappa*cexp(I*(gauge[hop[ix][nu]][nu-D] + phase));
         }

         mul_1_plus_gamma(nu,s+hop[ix][nu],&t);

         r[ix].s[0]+=c*t.s[0];
         r[ix].s[1]+=c*t.s[1];
      }
   }
}

