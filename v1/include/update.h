#ifndef UPDATE_H
#define UPDATE_H

#include "schwinger.h"

/* fermion.c */
extern double fermion_action1(spinor *phi,double kappa,double mu,double res);
extern double fermion_action2(spinor *phi,double kappa,double mu,double mu2,double res);
extern void fermion_force1(double f[V][D],spinor *phi,double kappa,double mu,double a,double res);
extern void fermion_force2(double f[V][D],spinor *phi,double kappa,double mu1,double mu2,double a,double res);
extern double setpf2(spinor *pf,double kappa,double mu1,double mu2,double res);
extern double setpf1(spinor *pf,double kappa,double mu1);

/* gauge.c   */
extern double gauge_action(double b);
extern void gauge_force(double f[V][D],double beta,double c);
extern double top_charge(void);

/* hmc.c */
extern double hmc(hmc_params_t *hpara,act_params_t *apara);
extern void reversibility(hmc_params_t *hpara,act_params_t *apara);

#endif
