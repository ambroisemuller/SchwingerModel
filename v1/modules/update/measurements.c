#include "measurements.h"
#include "update.h"
#include "linalg.h"
#include "dirac.h"
#include <stdlib.h>

double kappa, beta, eps, res_act, res_frc;
double *mus;
int nmax, npf;

char data_folder[128];
char plot_folder[128];
char study_folder[128];

FILE *file_dH;
const char *filename_dH_ = "dH.csv";
char filename_dH[64];

FILE *file_acc;
const char *filename_acc_ = "acc.csv";
char filename_acc[64];

double plaq;
FILE *file_plaq;
const char *filename_plaq_ = "plaq.csv";
char filename_plaq[64];

double qtop;
FILE *file_qtop;
const char *filename_qtop_ = "qtop.csv";
char filename_qtop[64];

double condensate;
FILE *file_condensate;
const char *filename_condensate_ = "condensate.csv";
char filename_condensate[64];

double corr_P__P_[V];   /*  < P(x) P(0) >       */
double corr_A0_P_[V];   /*  < A_0(x) P(0) >     */
double corr_A1_P_[V];   /*  < A_1(x) P(0) >     */
double corr_P__A0[V];   /*  < P(x) A_0(0) >     */
double corr_P__A1[V];   /*  < P(x) A_1(0) >     */
double corr_V0_V0[V];   /*  < V_0(x) V_0(0) >   */
double corr_V0_V1[V];   /*  < V_0(x) V_1(0) >   */
double corr_V1_V0[V];   /*  < V_1(x) V_0(0) >   */
double corr_V1_V1[V];   /*  < V_1(x) V_1(0) >   */

static spinor **phi;    /* workspace: CG solution spinor field */ 
static spinor *eta;     /* workspace: CG source spinor field   */ 
static spinor *temp;

FILE *file_PP;
const char *filename_PP_ = "PP.csv";
char filename_PP[64];
FILE *file_A0P;
const char *filename_A0P_ = "A0P.csv";
char filename_A0P[64];
FILE *file_A1P;
const char *filename_A1P_ = "A1P.csv";
char filename_A1P[64];
FILE *file_PA0;
const char *filename_PA0_ = "PA0.csv";
char filename_PA0[64];
FILE *file_PA1;
const char *filename_PA1_ = "PA1.csv";
char filename_PA1[64];
FILE *file_V0V0;
const char *filename_V0V0_ = "V0V0.csv";
char filename_V0V0[64];
FILE *file_V0V1;
const char *filename_V0V1_ = "V0V1.csv";
char filename_V0V1[64];
FILE *file_V1V0;
const char *filename_V1V0_ = "V1V0.csv";
char filename_V1V0[64];
FILE *file_V1V1;
const char *filename_V1V1_ = "V1V1.csv";
char filename_V1V1[64];


static void compute_2pt(){
    double mu = mus[0]; 
    int b_index;                    
    int pt1 = 0, pt2 = L*(T/2) + L/2;           
    complex double S00, S01, S10, S11; 
    complex double S00_, S01_, S10_, S11_; 

    /* compute S(x, y) */
    for (b_index=0; b_index<2; b_index++) {
        zero_spinor(phi[b_index]);                      /* Initialize solution to 0     */
        zero_spinor(eta);                               /* Initialize source to 0       */
        (eta[pt1]).s[b_index] = 1.;                     /* Set source to dirac delta    */
        gamma5(eta);                                    /* Apply gamma5 to source       */
        cg(eta, phi[b_index], kappa, mu, eps, nmax);    /* Solve the dirac equation     */
        dirac(phi[b_index], phi[b_index], kappa, mu);   /* Apply Dirac operator         */
        gamma5(phi[b_index]);                           /* Apply gamma5 to solution     */
    }
    S00 = phi[0][pt2].s[0];
    S01 = phi[1][pt2].s[0];
    S10 = phi[0][pt2].s[1];
    S11 = phi[1][pt2].s[1];

    /* compute S(y, x) */
    for (b_index=0; b_index<2; b_index++) {
        zero_spinor(phi[b_index]);                      /* Initialize solution to 0     */
        zero_spinor(eta);                               /* Initialize source to 0       */
        (eta[pt2]).s[b_index] = 1.;                     /* Set source to dirac delta    */
        gamma5(eta);                                    /* Apply gamma5 to source       */
        cg(eta, phi[b_index], kappa, mu, eps, nmax);    /* Solve the dirac equation     */
        dirac(phi[b_index], phi[b_index], kappa, mu);   /* Apply Dirac operator         */
        gamma5(phi[b_index]);                           /* Apply gamma5 to solution     */
    }
    S00_ = phi[0][pt1].s[0];
    S01_ = phi[1][pt1].s[0];
    S10_ = phi[0][pt1].s[1];
    S11_ = phi[1][pt1].s[1];

    printf("\nS(middle, zero)\n%4.3e + i%4.3e   %4.3e + i%4.3e\n%4.3e + i%4.3e   %4.3e + i%4.3e\n", creal(S00), cimag(S00), creal(S01), cimag(S01), creal(S10), cimag(S10), creal(S11), cimag(S11));
    printf("\nS(zero, middle)\n%4.3e + i%4.3e   %4.3e + i%4.3e\n%4.3e + i%4.3e   %4.3e + i%4.3e\n", creal(S00_), cimag(S00_), creal(S01_), cimag(S01_), creal(S10_), cimag(S10_), creal(S11_), cimag(S11_));

}


static void compute_correlators()
{
    /* initialize */
    double mu = mus[0]; /* or just 0? */
    int b_index, site_index;                            /* Spin and lattice indices     */
    complex double S00_S11c, S01_S10c, S00_S01c, S10_S11c, S00_S10c, S01_S11c;
    int zero_index, x_index, t_index;                   /* zero = 0, otherwise middle of lattice     L*(T/2) + L/2       */
    int x, t, actual_x_index, actual_t_index, actual_index;           /* workspace */
    int sample_index;
    
    /* zero all global variables */
    condensate = 0;
    
    for (site_index = 0; site_index < V; site_index++) {
        corr_P__P_[site_index] = 0;
        corr_A0_P_[site_index] = 0; 
        corr_A1_P_[site_index] = 0;
        corr_P__A0[site_index] = 0;
        corr_P__A1[site_index] = 0;
        corr_V0_V0[site_index] = 0;
        corr_V0_V1[site_index] = 0;
        corr_V1_V0[site_index] = 0;
        corr_V1_V1[site_index] = 0;
    }

    for (sample_index=0; sample_index<CORR_SAMPLES; sample_index++)
    {
        t_index = rand()%T;
        x_index = rand()%L;
        zero_index = t_index*L + x_index;

        printf("(%d, %d)", t_index, x_index);

        /* compute S(x, 0) */
        for (b_index=0; b_index<2; b_index++) {
            zero_spinor(phi[b_index]);                      /* Initialize solution to 0     */
            zero_spinor(eta);                               /* Initialize source to 0       */
            zero_spinor(temp);
            (eta[zero_index]).s[b_index] = 1.;              /* Set source to dirac delta    */
            gamma5(eta);                                    /* Apply gamma5 to source       */
            cg(eta, temp, kappa, mu, eps, nmax);            /* Solve the dirac equation     */
            dirac(temp, phi[b_index], kappa, mu);           /* Apply Dirac operator         */
            gamma5(phi[b_index]);                           /* Apply gamma5 to solution     */
        };

        /* fermion condensate */
        condensate = creal(phi[0][zero_index].s[0] + phi[1][site_index].s[1]); 

        /* measure correlators */
        for (site_index = 0; site_index < V; site_index++) 
        {
            
            t = site_index / L;
            x = site_index % L;
            actual_t_index = (t-t_index+T)%T;
            actual_x_index = (x-x_index+L)%L;
            actual_index = actual_t_index*L + actual_x_index;

            corr_P__P_[actual_index] += creal(
                phi[0][site_index].s[0] * conj(phi[0][site_index].s[0])
                + phi[0][site_index].s[1] * conj(phi[0][site_index].s[1])
                + phi[1][site_index].s[0] * conj(phi[1][site_index].s[0])
                + phi[1][site_index].s[1] * conj(phi[1][site_index].s[1])
            )/((double)(CORR_SAMPLES));
        
            S00_S01c = phi[0][site_index].s[0] * conj(phi[1][site_index].s[0]);
            S10_S11c = phi[0][site_index].s[1] * conj(phi[1][site_index].s[1]);
            S00_S10c = phi[0][site_index].s[0] * conj(phi[0][site_index].s[1]);
            S01_S11c = phi[1][site_index].s[0] * conj(phi[1][site_index].s[1]);

            corr_A0_P_[actual_index] += (2*creal(S00_S01c + S10_S11c))/((double)(CORR_SAMPLES)); 
            corr_A1_P_[actual_index] += (2*cimag(S00_S01c + S10_S11c))/((double)(CORR_SAMPLES));
            corr_P__A0[actual_index] += (- 2*creal(S00_S10c + S01_S11c))/((double)(CORR_SAMPLES));
            corr_P__A1[actual_index] += (2*cimag(S00_S10c + S01_S11c))/((double)(CORR_SAMPLES));

            S00_S11c = phi[0][site_index].s[0] * conj(phi[1][site_index].s[1]);
            S01_S10c = phi[1][site_index].s[0] * (phi[0][site_index].s[1]);

            corr_V0_V0[actual_index] += /* creal(S00_S11c); */ (2*creal(S00_S11c) - 2*creal(S01_S10c))/((double)(CORR_SAMPLES)); /* something must be wrong here, no consistent with G00 = -G11 result */
            corr_V0_V1[actual_index] += /* cimag(S00_S11c); */ (- 2*cimag(S00_S11c) + 2*cimag(S01_S10c))/((double)(CORR_SAMPLES));
            corr_V1_V0[actual_index] += /* cimag(S01_S10c); */ (+ 2*cimag(S00_S11c) + 2*cimag(S01_S10c))/((double)(CORR_SAMPLES));
            corr_V1_V1[actual_index] += /* creal(S01_S10c); */ (+ 2*creal(S00_S11c) + 2*creal(S01_S10c))/((double)(CORR_SAMPLES));
            
        }
    }
    printf("\n");

}

void datafile_headers(hmc_params_t *hmc_params,act_params_t *act_params)
{
    int i;
    int b_index;
    char bc[2];
    char command[256];

    kappa = act_params->kappa;
    beta = act_params->beta;
    npf = hmc_params->npf;
    mus = hmc_params->mu;
    res_act = hmc_params->res_act;
    res_frc = hmc_params->res_frc;
    eps = res_frc;
    nmax = 1e3;

    if (PERIODIC == 0){
        bc[0] = 'A';
    } else {
        bc[0] = 'P';
    }
    bc[1] = '\0';
    
    sprintf(data_folder, "%sT%d_L%d_%s", DATA_FOLDER, T, L, bc);
    sprintf(command, "mkdir %s", data_folder);
    system(command);
    sprintf(data_folder, "%sT%d_L%d_%s/k%.4f_b%.3f/", DATA_FOLDER, T, L, bc, kappa, beta);
    sprintf(study_folder, "T%d_L%d_%s/k%.4f_b%.3f/", T, L, bc, kappa, beta);
    sprintf(command, "mkdir %s", data_folder);
    system(command);
    sprintf(command, "cp ./infile_qed %s", data_folder);
    system(command);
    sprintf(data_folder, "%sT%d_L%d_%s/k%.4f_b%.3f/%s/", DATA_FOLDER, T, L, bc, kappa, beta, "observables");
    sprintf(plot_folder, "%sT%d_L%d_%s/k%.4f_b%.3f/%s/", DATA_FOLDER, T, L, bc, kappa, beta, "plots");
    sprintf(command, "mkdir %s", data_folder);
    system(command);
    sprintf(command, "mkdir %s", plot_folder);
    system(command);

    printf(data_folder);
    printf(plot_folder);

    sprintf(filename_dH, "%s%s", data_folder, filename_dH_);
    file_dH = fopen(filename_dH, "w");
    fprintf(file_dH, "time,dH\n");
    fclose(file_dH);

    sprintf(filename_acc, "%s%s", data_folder, filename_acc_);
    file_acc = fopen(filename_acc, "w");
    fprintf(file_acc, "time,acc\n");
    fclose(file_acc);

    sprintf(filename_plaq, "%s%s", data_folder, filename_plaq_);
    file_plaq = fopen(filename_plaq, "w");
    fprintf(file_plaq, "time,plaq\n");
    fclose(file_plaq);

    sprintf(filename_qtop, "%s%s", data_folder, filename_qtop_);
    file_qtop = fopen(filename_qtop, "w");
    fprintf(file_qtop, "time,qtop\n");
    fclose(file_qtop);

    sprintf(filename_condensate, "%s%s", data_folder, filename_condensate_);
    file_condensate = fopen(filename_condensate, "w");
    fprintf(file_condensate, "time,condensate\n");
    fclose(file_condensate);

    /* allocate workspace */
    phi = malloc(2*sizeof(spinor*));
    for (b_index=0; b_index<2; b_index++) {
        phi[b_index] = malloc(V*sizeof(spinor));
    }
    eta = malloc(V*sizeof(spinor));
    temp = malloc(V*sizeof(spinor));
    /* PP */
    sprintf(filename_PP, "%s%s", data_folder, filename_PP_);
    file_PP = fopen(filename_PP, "w");
    fprintf(file_PP, "time,");
    for (i=0; i<V; i++) { 
        fprintf(file_PP, "%d", i);
        if (i != V-1) { fprintf(file_PP, ",");}
    }
    fprintf(file_PP, "\n");
    fclose(file_PP);
    /* A0P */
    sprintf(filename_A0P, "%s%s", data_folder, filename_A0P_);
    file_A0P = fopen(filename_A0P, "w");
    fprintf(file_A0P, "time,");
    for (i=0; i<V; i++) { 
        fprintf(file_A0P, "%d", i);
        if (i != V-1) { fprintf(file_A0P, ",");}
    }
    fprintf(file_A0P, "\n");
    fclose(file_A0P);
    /* A1P */
    sprintf(filename_A1P, "%s%s", data_folder, filename_A1P_);
    file_A1P = fopen(filename_A1P, "w");
    fprintf(file_A1P, "time,");
    for (i=0; i<V; i++) { 
        fprintf(file_A1P, "%d", i);
        if (i != V-1) { fprintf(file_A1P, ",");}
    }
    fprintf(file_A1P, "\n");
    fclose(file_A1P);
    /* PA0 */
    sprintf(filename_PA0, "%s%s", data_folder, filename_PA0_);
    file_PA0 = fopen(filename_PA0, "w");
    fprintf(file_PA0, "time,");
    for (i=0; i<V; i++) { 
        fprintf(file_PA0, "%d", i);
        if (i != V-1) { fprintf(file_PA0, ",");}
    }
    fprintf(file_PA0, "\n");
    fclose(file_PA0);
    /* PA1 */
    sprintf(filename_PA1, "%s%s", data_folder, filename_PA1_);
    file_PA1 = fopen(filename_PA1, "w");
    fprintf(file_PA1, "time,");
    for (i=0; i<V; i++) { 
        fprintf(file_PA1, "%d", i);
        if (i != V-1) { fprintf(file_PA1, ",");}
    }
    fprintf(file_PA1, "\n");
    fclose(file_PA1);
    /* V0V0 */
    sprintf(filename_V0V0, "%s%s", data_folder, filename_V0V0_);
    file_V0V0 = fopen(filename_V0V0, "w");
    fprintf(file_V0V0, "time,");
    for (i=0; i<V; i++) { 
        fprintf(file_V0V0, "%d", i);
        if (i != V-1) { fprintf(file_V0V0, ",");}
    }
    fprintf(file_V0V0, "\n");
    fclose(file_V0V0);
    /* V0V1 */
    sprintf(filename_V0V1, "%s%s", data_folder, filename_V0V1_);
    file_V0V1 = fopen(filename_V0V1, "w");
    fprintf(file_V0V1, "time,");
    for (i=0; i<V; i++) { 
        fprintf(file_V0V1, "%d", i);
        if (i != V-1) { fprintf(file_V0V1, ",");}
    }
    fprintf(file_V0V1, "\n");
    fclose(file_V0V1);
    /* V1V0 */
    sprintf(filename_V1V0, "%s%s", data_folder, filename_V1V0_);
    file_V1V0 = fopen(filename_V1V0, "w");
    fprintf(file_V1V0, "time,");
    for (i=0; i<V; i++) { 
        fprintf(file_V1V0, "%d", i);
        if (i != V-1) { fprintf(file_V1V0, ",");}
    }
    fprintf(file_V1V0, "\n");
    fclose(file_V1V0);
    /* V1V1 */
    sprintf(filename_V1V1, "%s%s", data_folder, filename_V1V1_);
    file_V1V1 = fopen(filename_V1V1, "w");
    fprintf(file_V1V1, "time,");
    for (i=0; i<V; i++) { 
        fprintf(file_V1V1, "%d", i);
        if (i != V-1) { fprintf(file_V1V1, ",");}
    }
    fprintf(file_V1V1, "\n");
    fclose(file_V1V1);

}

void measure_and_record(spinor** pf, double time, double dH, double acc)
{

    int i;

    file_dH = fopen(filename_dH, "a");
    fprintf(file_dH, "%4.3e,%4.3e\n", time, dH);
    fclose(file_dH);

    file_acc = fopen(filename_acc, "a");
    fprintf(file_acc, "%4.3e,%4.3e\n", time, acc);
    fclose(file_acc);

    plaq = gauge_action(-1)/V;
    file_plaq = fopen(filename_plaq, "a");
    fprintf(file_plaq, "%4.3e,%4.3e\n", time, plaq);
    fclose(file_plaq);

    qtop = top_charge();
    file_qtop = fopen(filename_qtop, "a");
    fprintf(file_qtop, "%4.3e,%4.3e\n", time, qtop);
    fclose(file_qtop);

    compute_correlators();

    file_condensate = fopen(filename_condensate, "a");
    fprintf(file_condensate, "%4.3e,%4.3e\n", time, condensate);
    fclose(file_condensate);

    /* PP */
    file_PP = fopen(filename_PP, "a");
    fprintf(file_PP, "%4.3e,", time);
    for (i=0; i<V; i++) { 
        fprintf(file_PP, "%4.3e", corr_P__P_[i]); 
        if (i != V-1) { fprintf(file_PP, ",");}
    }
    fprintf(file_PP, "\n");
    fclose(file_PP);
    /* A0P */
    file_A0P = fopen(filename_A0P, "a");
    fprintf(file_A0P, "%4.3e,", time);
    for (i=0; i<V; i++) { 
        fprintf(file_A0P, "%4.3e", corr_A0_P_[i]); 
        if (i != V-1) { fprintf(file_A0P, ",");}
    }
    fprintf(file_A0P, "\n");
    fclose(file_A0P);
    /* A1P */
    file_A1P = fopen(filename_A1P, "a");
    fprintf(file_A1P, "%4.3e,", time);
    for (i=0; i<V; i++) { 
        fprintf(file_A1P, "%4.3e", corr_A1_P_[i]); 
        if (i != V-1) { fprintf(file_A1P, ",");}
    }
    fprintf(file_A1P, "\n");
    fclose(file_A1P);
    /* PA0 */
    file_PA0 = fopen(filename_PA0, "a");
    fprintf(file_PA0, "%4.3e,", time);
    for (i=0; i<V; i++) { 
        fprintf(file_PA0, "%4.3e", corr_P__A0[i]); 
        if (i != V-1) { fprintf(file_PA0, ",");}
    }
    fprintf(file_PA0, "\n");
    fclose(file_PA0);
    /* PA1 */
    file_PA1 = fopen(filename_PA1, "a");
    fprintf(file_PA1, "%4.3e,", time);
    for (i=0; i<V; i++) { 
        fprintf(file_PA1, "%4.3e", corr_P__A1[i]); 
        if (i != V-1) { fprintf(file_PA1, ",");}
    }
    fprintf(file_PA1, "\n");
    fclose(file_PA1);
    /* V0V0 */
    file_V0V0 = fopen(filename_V0V0, "a");
    fprintf(file_V0V0, "%4.3e,", time);
    for (i=0; i<V; i++) { 
        fprintf(file_V0V0, "%4.3e", corr_V0_V0[i]); 
        if (i != V-1) { fprintf(file_V0V0, ",");}
    }
    fprintf(file_V0V0, "\n");
    fclose(file_V0V0);
    /* V0V1 */
    file_V0V1 = fopen(filename_V0V1, "a");
    fprintf(file_V0V1, "%4.3e,", time);
    for (i=0; i<V; i++) { 
        fprintf(file_V0V1, "%4.3e", corr_V0_V1[i]); 
        if (i != V-1) { fprintf(file_V0V1, ",");}
    }
    fprintf(file_V0V1, "\n");
    fclose(file_V0V1);
    /* V1V0 */
    file_V1V0 = fopen(filename_V1V0, "a");
    fprintf(file_V1V0, "%4.3e,", time);
    for (i=0; i<V; i++) { 
        fprintf(file_V1V0, "%4.3e", corr_V1_V0[i]); 
        if (i != V-1) { fprintf(file_V1V0, ",");}
    }
    fprintf(file_V1V0, "\n");
    fclose(file_V1V0);
    /* V1V1 */
    file_V1V1 = fopen(filename_V1V1, "a");
    fprintf(file_V1V1, "%4.3e,", time);
    for (i=0; i<V; i++) { 
        fprintf(file_V1V1, "%4.3e", corr_V1_V1[i]); 
        if (i != V-1) { fprintf(file_V1V1, ",");}
    }
    fprintf(file_V1V1, "\n");
    fclose(file_V1V1);
}

void run_plot_scripts()
{
    char command[256];
    int out;
    
    sprintf(command, "cd ../studies/ && python3 plot.py %s", study_folder);
    out = system("");
    out = system(command);
    if (out == 0){
        printf("plotting of general observables successful\n");
    }
    sprintf(command, "cd ../studies/ && python3 plot_corr.py %s", study_folder);
    out = system("");
    out = system(command);
    if (out == 0){
        printf("plotting of correlators successful\n");
    }
    sprintf(command, "cd ../studies/ && python3 pcac_mass.py %s", study_folder);
    out = system("");
    out = system(command);
    if (out == 0){
        printf("plotting of PCAC mass successful\n");
    }
    sprintf(command, "cd ../studies/ && python3 pion_mass.py %s", study_folder);
    out = system("");
    out = system(command);
    if (out == 0){
        printf("plotting of pion mass successful\n");
    }
}