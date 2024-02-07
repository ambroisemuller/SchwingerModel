/*********************************************************
 * 
 *  File hopping.c
 *
 *  Initialization of the hopping field for a D dimensional
 *  lattice of size V=L**D
 *  The index of a point with coordinates (n_0,n_1,..,n_{d-1})
 *  is i=sum_k n_k L**k
 *  The index of its neighbor in the positive direction nu 
 *  is hop[i][mu]
 *  In the negative direction it is hop[i][D+mu]
 *  Periodic boundary conditions are imposed on the points
 *
 ********************************************************/

#include "geometry.h"

/**
 * For a rectangular lattice, but only in 2D
*/
void hopping(int hop[V][2*D]) 
{
   int i, x, t;
   for (i=0; i<V; i++){
      x = i % L;
      t = i / L;
      hop[i][0] = t*L + ((x+1) % L);
      hop[i][1] = ((t+1) % T)*L + x;
      hop[i][2] = t*L + ((x-1+L) % L);
      hop[i][3] = ((t-1+T) % T)*L + x;

   }
}

#ifdef SQUARE_LATTICE

/**
 * For a square lattice only, but in arbitrary dimension
*/
void hopping(int hop[V][2*D])
{
   int x,y,Lk;
   int xk,k,dxk;

   /* go through all the points*/
   
   for (x=0;x<V;x++)
   {
      Lk=V;
      y =x;

      /* go through the components k*/
      for (k=D-1;k>=0;k--)
      {
         Lk/=L;                        /* pow(L,k)      */
         xk =y/Lk;                     /* kth component */
         y  =y-xk*Lk;                  /* y<-y%Lk       */

         /* forward */
         if (xk<L-1) 
            dxk=Lk;
         else        
            dxk=Lk*(1-L);
         hop[x][k]=x+dxk;

         /* backward */
         if (xk>0)   
            dxk=-Lk;
         else        
            dxk=Lk*(L-1);
         hop[x][k+D]=x+dxk;
      }
   }
} 

#endif