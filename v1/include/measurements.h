#ifndef MEASUREMENTS_H
#define MEASUREMENTS_H

#include "schwinger.h"

/* Standard */
#define MEASURE_DH      1
#define MEASURE_ACC     1
#define MEASURE_PLAQ    1
#define MEASURE_QTOP    1

/*  Correlators */
#define MEASURE_CORR    1

#define DATA_FOLDER     "../studies/"

/* measurements.c */
extern void datafile_headers(hmc_params_t *hmc_params,act_params_t *act_params);
extern void measure_and_record(spinor** pf, double time, double dH, int acc);
extern void run_plot_scripts(void);

#endif