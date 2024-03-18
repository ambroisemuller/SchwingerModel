#ifndef MEASUREMENTS_H
#define MEASUREMENTS_H

#include "schwinger.h"

#define DATA_FOLDER     "../studies/"

/* measurements.c */
extern void datafile_headers(hmc_params_t *hmc_params,act_params_t *act_params);
extern void measure_and_record(spinor** pf, double time, double dH, double acc);
extern void run_plot_scripts(void);

#endif