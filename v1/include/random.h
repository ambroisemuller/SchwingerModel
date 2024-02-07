#ifndef RANDOM_H
#define RANDOM_H

#include "schwinger.h"

/*gauss_rand.c*/
extern void gauss_rand(int n,double *rand);
extern void gauss_spinor(spinor *rand);

/* ranlxd.c */
extern void ranlxd(double r[],int n);
extern void rlxd_init(int level,int seed);
extern int rlxd_size(void);
extern void rlxd_get(int state[]);
extern void rlxd_reset(int state[]);

#endif
