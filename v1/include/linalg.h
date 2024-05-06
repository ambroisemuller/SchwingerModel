#ifndef LINALG_H
#define LINALG_H

#include "schwinger.h"

/* cg.c */
extern void cg(spinor *s,spinor *r,double mass,double mu,double eps,int nmax);

/* linalg.c */
extern void zero_frc(double frc[V][D]);
extern void add_frc(double frc[V][D],double frc2[V][D]);
extern double norm_frc(double frc[V][D]);
extern void assign_link(double s[V][D],double r[V][D]);
extern double cmp_link(double s[V][D],double r[V][D]);
extern void assign_spinor(spinor *s,spinor *r);
extern void zero_spinor(spinor *r);
extern double complex scalar_prod(spinor *s,spinor *r);
extern double complex local_scalar_prod(spinor phi,spinor psi);
extern void mulr_add_assign(spinor *c,spinor *a,double z,spinor *b);
extern void mulc_add_assign(spinor *a,spinor *b,double complex z);

#define MAX_CG_STEPS 10000

#endif
