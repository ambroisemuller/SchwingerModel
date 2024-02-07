/* 
 *    File gauss_rand.c
 *
 *    void gauss_rand(int n,double *rand)
 *        Generates n double precision gaussian random numbers
 *        Box-Muller procedure
 *
 *    void gauss_spinor(spinor *s)
 *        Generates a gaussian distributed spinor field over 
 *        the entire lattice volume
 */

#include "random.h"


void gauss_rand(int n,double *rand)
{
   static int called=0;
   static double tpi;
   int i;
   double tmp[2],r,phi;

   if (!called) 
   {
      tpi=8.*atan(1.);
      called=1;
   }

   i=0;
   while(i<n)
   {
      /* two random numbers, flat distribution in [0,1) */
      ranlxd(tmp,2);

      /* compute polar coordinates: angle and radius */
      phi=tmp[0]*tpi;
      r  =sqrt(-log(1.-tmp[1])); /* map second number [0,1) -> (0,1] */

      rand[i]=r*cos(phi); 
      i++;
      /*compute second only if requested */
      if (i<n) rand[i]=r*sin(phi); 
      i++;
   }
}


void gauss_spinor(spinor *s)
{
   int ix;
   double rv[4];

   for (ix=0;ix<V;ix++)
   {
      gauss_rand(4,rv);
      s[ix].s[0]=rv[0]+I*rv[1];
      s[ix].s[1]=rv[2]+I*rv[3];
   }
}

