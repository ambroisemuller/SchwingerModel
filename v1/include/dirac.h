#ifndef DIRAC_H
#define DIRAC_H

#include "schwinger.h"

/* dirac.c */
extern void gamma5(spinor *s);
extern void mul_1_plus_gamma(int mu,spinor *s,spinor *r);
extern void dirac(spinor *s,spinor *r,double mass,double mu);

#endif
