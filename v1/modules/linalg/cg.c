/*******************************************************************************
 *
 * File cg.c
 *
 * Copyright (C) 2017 Stefan Schaefer
 *
 * This software is distributed under the terms of the GNU General Public
 * License (GPL)
 *
 * Conjugate gradient to solve the Dirac equation
 *
 * The externally accessible functions are
 *
 * void cg(spinor *s,spinor *r,double kappa,double mu,double eps,int nmax)
 *   Conjugate gradient algorithm to find solution to 
 *   
 *   (D^dag*D + mu^2) r = s
 *
 *   kappa  hopping parameter in fermion matrix M
 *   eps    norm of the residue
 *   nmax   maximal number of iterations
 *
 *******************************************************************************/

#include "linalg.h"
#include "dirac.h"

static spinor work[V],p[V],res[V],ap[V];


static void mxv(spinor *s,spinor *r,double kappa,double mu)
{
   dirac(s,work,kappa,mu);
   gamma5(work);
   dirac(work,r,kappa,-mu);
   gamma5(r);
}


void cg(spinor *s,spinor *r,double kappa,double mu,double eps,int nmax)
{
   double rsold,rsnew,alpha;
   double complex c;
   int i;

   zero_spinor(r);
   assign_spinor(s,res);
   assign_spinor(res,p);
   rsold=creal(scalar_prod(res,res));

   for (i=0;i<nmax;i++)
   {
      mxv(p,ap,kappa,mu);
      c=scalar_prod(p,ap);
      alpha=rsold/creal(c);
      mulr_add_assign(r,r,alpha,p);
      mulr_add_assign(res,res,-alpha,ap);
      rsnew=creal(scalar_prod(res,res));
      if (sqrt(rsnew)<eps)
         break;
      mulr_add_assign(p,res,rsnew/rsold,p);
      rsold=rsnew;
   }

   if (i==nmax) 
   {
      fprintf(stderr,"CG did not converge");
      exit(1);
   }
}

