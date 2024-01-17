/*******************************************************************************
 *
 * File linalg.c
 *
 * Copyright (C) 2017 Stefan Schaefer
 *
 * This software is distributed under the terms of the GNU General Public
 * License (GPL)
 *
 * Helper routines for spinor and forces
 *
 * The externally accessible functions are
 *
 * void zero_frc(double frc[V][D])
 *   frc=0;
 *
 * void add_frc(double a[V][D],double b[V][D])
 *   a+=b
 *
 * double norm_frc(double frc[V][D])
 *   returns |frc|^2
 *
 * void assign_link(double s[V][D],double r[V][D])
 *   r=s
 *
 * double cmp_link(double s[V][D],double r[V][D])
 *   Returns the max over all lattice points and directions of the absolute
 *   deviation between the elements of r and s
 *
 * void assign_spinor(spinor *s,spinor *r)
 *   r=s
 *
 * void zero_spinor(spinor *r)
 *   r=0
 *
 * double complex scalar_prod(spinor *s,spinor *r)
 *   returns (s,r)
 * 
 * void mulr_add_assign(spinor *c,spinor *a,double z,spinor *b)
 *   c=a+z*b
 *
 * void mulc_add_assign(spinor *a,spinor *b,double complex z)
 *   a=a+z*b
 *
 *******************************************************************************/

#include "linalg.h"


void zero_frc(double frc[V][D])
{
   int i,mu;

   for (i=0;i<V;i++) 
      for (mu=0;mu<D;mu++) 
         frc[i][mu]=0.0;
}


void add_frc(double a[V][D],double b[V][D])
{
   int i,mu;

   for (i=0;i<V;i++) 
      for (mu=0;mu<D;mu++) 
         a[i][mu]+=b[i][mu];
}


double norm_frc(double frc[V][D])
{
   int i,mu;
   double norm=0;

   for (i=0;i<V;i++) 
      for (mu=0;mu<D;mu++) 
         norm+=frc[i][mu]*frc[i][mu];

   return norm;
}


void assign_link(double s[V][D],double r[V][D])
{
   int i,mu;

   for (i=0;i<V;i++) 
      for (mu=0;mu<D;mu++)
         r[i][mu]=s[i][mu];
}


double cmp_link(double s[V][D],double r[V][D])
{
   int i,mu;
   double dev,dmax;

   dmax=0.0;

   for(i=0;i<V;i++)
   {   
      for(mu=0;mu<D;mu++)
      {
         dev=fabs(r[i][mu]-s[i][mu]);
         if (dev>dmax)
            dmax=dev;
      }
   }

   return dmax;
}   


void assign_spinor(spinor *s,spinor *r)
{
   int i;

   for (i=0;i<V;i++)
      r[i]=s[i];
}


void zero_spinor(spinor *r)
{
   int i;

   for (i=0;i<V;i++)
   {
      r[i].s[0]=0;
      r[i].s[1]=0;
   }
}


double complex scalar_prod(spinor *s,spinor *r)
{
   int i;
   double complex x=0.;

   for (i=0;i<V;i++)
   {
      x+=conj(s[i].s[0])*r[i].s[0];
      x+=conj(s[i].s[1])*r[i].s[1];
   }

   return x;
}

double complex local_scalar_prod(spinor phi,spinor psi)
{
   double complex x=0.;

   x+=conj(phi.s[0])*psi.s[0];
   x+=conj(phi.s[1])*psi.s[1];

   return x;
}


void mulr_add_assign(spinor *c,spinor *a,double z,spinor *b)
{
   int i;

   for (i=0;i<V;i++)
   {
      c[i].s[0]=a[i].s[0]+z*b[i].s[0];
      c[i].s[1]=a[i].s[1]+z*b[i].s[1];
   }
}


void mulc_add_assign(spinor *a,spinor *b,double complex z)
{
   int i;

   for (i=0;i<V;i++)
   {
      a[i].s[0]+=z*b[i].s[0];
      a[i].s[1]+=z*b[i].s[1];
   }
}

