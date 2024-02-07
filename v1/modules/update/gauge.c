/*******************************************************************************
 *
 * File gauge.c
 *
 * Copyright (C) 2017 Stefan Schaefer, 2021 Mattia Dalla Brida
 *
 * This software is distributed under the terms of the GNU General Public
 * License (GPL)
 *
 * For the 2D Schwinger model: gauge action, gauge force, and topological
 * charge
 *
 * The externally accessible functions are
 *
 * double gauge_action(double beta)
 *   Returns the gauge action for coupling parameter beta
 *   S=-beta*sum_x real(U(x,0)U(x+0,1)U^*(x+1,0)U^*(x,1))
 *
 * void gauge_force(double f[V][D],double beta,double a)
 *   Sets the force field f to a times the derivative of the
 *   gauge action with parameter beta
 *
 * double top_charge(void)
 *   Returns the geometricial definition of the topological charge:
 *   Q=(1/2*PI)*sum_x Im(ln(U_01)), where U_01 is the plaquette in
 *   the 0-1 plane
 *
 *******************************************************************************/

#include "update.h"


double gauge_action(double beta)
{
   int i,imu,inu;
   double act=0.0,phi;

   for (i=0;i<V;i++)
   {
      imu=hop[i][0];
      inu=hop[i][1];
      phi=gauge[i][0]+gauge[imu][1]-gauge[inu][0]-gauge[i][1];
      act+=cos(phi);
   }

   return -beta*act;
}


void gauge_force(double f[V][D],double beta,double a)
{
   int ix,mu,nu,ix1,ix2;

   for (ix=0;ix<V;ix++)
   {
      for (mu=0;mu<D;mu++)
      {
         for (nu=0;nu<D;nu++)
         {
            if (mu==nu) continue;
            /*forward*/
            ix1=hop[ix][mu];
            ix2=hop[ix][nu];
            f[ix][mu]=-a*beta*sin(gauge[ix][mu]+gauge[ix1][nu]-gauge[ix2][mu]-gauge[ix][nu]);
            /*backward*/
            ix1=hop[ix][D+nu];
            ix2=hop[ix1][mu];
            f[ix][mu]+=-a*beta*sin(gauge[ix][mu]-gauge[ix2][nu]-gauge[ix1][mu]+gauge[ix1][nu]);
         }
      }
   }
}


double top_charge(void)
{
   static int called=0;
   static double tpi;
   int i,imu,inu;
   double qtop=0.0,offset,phi;

   if (!called)
   {
      tpi=8.*atan(1.);
      called=1;
   }

   for (i=0;i<V;i++)
   {
      imu=hop[i][0];
      inu=hop[i][1];
      phi=gauge[i][0]+gauge[imu][1]-gauge[inu][0]-gauge[i][1];
      offset=floor((phi+0.5*tpi)/tpi)*tpi;
      qtop+=phi-offset;
   }

  return qtop/tpi;
}

