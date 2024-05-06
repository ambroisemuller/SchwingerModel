#ifndef LATTICE_H
#define LATTICE_H

/* Dimension of the lattice */
#define D 2
/* spatial extent of the lattice */
#define L 128
/* temporal extent of the lattice */
#define T 256
/* lattice volume, needs to be adjusted according to number of dimensions*/
#define V (L*T)

/* Periodicity of pseudofermion fields in time */
#define PERIODIC 0

#ifdef CONTROL 
#define EXTERN 
#undef CONTROL
#else
#define EXTERN extern
#endif

EXTERN double gauge[V][D];
EXTERN int hop[V][2*D];

#endif
