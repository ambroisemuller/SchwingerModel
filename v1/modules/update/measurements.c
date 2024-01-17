#include "measurements.h"
#include "update.h"
#include "linalg.h"
#include "dirac.h"

double kappa, beta, eps, res_act, res_frc;
double *mus;
int nmax, npf;

#if MEASURE_DH
    FILE *file_dH;
    const char *filename_dH_ = "dH.csv";
    char filename_dH[64];
    #endif
#if MEASURE_ACC
    FILE *file_acc;
    const char *filename_acc_ = "acc.csv";
    char filename_acc[64];
#endif
#if MEASURE_PLAQ
    double plaq;
    FILE *file_plaq;
    const char *filename_plaq_ = "plaq.csv";
    char filename_plaq[64];
#endif
#if MEASURE_QTOP
    double qtop;
    FILE *file_qtop;
    const char *filename_qtop_ = "qtop.csv";
    char filename_qtop[64];
#endif

#if MEASURE_CORR
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
    FILE *file_PP;
    const char *filename_PP_ = "PP.csv";
    char filename_PP[64];
    FILE *file_A0P;
    const char *filename_A0P_ = "A0P.csv";
    char filename_A0P[64];
    FILE *file_A1P;
    const char *filename_A1P_ = "A1P.csv";
    char filename_A1P[64];
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
    double mu = mus[0]; 
    int b_index, site_index;                            /* Spin and lattice indices     */
    int zero_index = L*(T/2) + L/2;                     /* middle of lattice            */
    complex double S00_S11c, S01_S10c, S00_S01c, S10_S11c, S00_S10c, S01_S11c;

    /* compute S(x, 0) */
    for (b_index=0; b_index<2; b_index++) {
        zero_spinor(phi[b_index]);                      /* Initialize solution to 0     */
        zero_spinor(eta);                               /* Initialize source to 0       */
        (eta[zero_index]).s[b_index] = 1.;              /* Set source to dirac delta    */
        gamma5(eta);                                    /* Apply gamma5 to source       */
        cg(eta, phi[b_index], kappa, mu, eps, nmax);    /* Solve the dirac equation     */
        dirac(phi[b_index], phi[b_index], kappa, mu);   /* Apply Dirac operator         */
        gamma5(phi[b_index]);                           /* Apply gamma5 to solution     */
    }

    /* measure correlators */
    for (site_index = 0; site_index < V; site_index++) {
        corr_P__P_[site_index] = creal(
            phi[0][site_index].s[0] * conj(phi[0][site_index].s[0])
            + phi[0][site_index].s[1] * conj(phi[0][site_index].s[1])
            + phi[1][site_index].s[0] * conj(phi[1][site_index].s[0])
            + phi[1][site_index].s[1] * conj(phi[1][site_index].s[1])
        );
        S00_S01c = phi[0][site_index].s[0] * conj(phi[1][site_index].s[0]);
        S10_S11c = phi[0][site_index].s[1] * conj(phi[1][site_index].s[1]);
        S00_S10c = phi[0][site_index].s[0] * conj(phi[0][site_index].s[1]);
        S01_S11c = phi[1][site_index].s[0] * conj(phi[1][site_index].s[1]);
        corr_A0_P_[site_index] = 2*creal(- S00_S10c - S01_S11c); /* S00_S01c + S10_S11c ; */
        corr_A1_P_[site_index] = 2*cimag(- S00_S10c - S01_S11c); /* S00_S01c + S10_S11c ; */
        /*
        corr_P__A0[site_index] = 
        );
        corr_P__A1[site_index] = 
        );
        */
        S00_S11c = phi[0][site_index].s[0] * conj(phi[1][site_index].s[1]);
        S01_S10c = phi[1][site_index].s[0] * conj(phi[0][site_index].s[1]);
        corr_V0_V0[site_index] = 2*creal(S00_S11c) - 2*creal(S01_S10c);
        corr_V0_V1[site_index] = - 2*cimag(S00_S11c) + 2*cimag(S01_S10c);
        corr_V1_V0[site_index] = 2*cimag(S00_S11c) + 2*cimag(S01_S10c);
        corr_V1_V1[site_index] = 2*creal(S00_S11c) + 2*creal(S01_S10c);
    }
}
#endif

void datafile_headers(hmc_params_t *hmc_params,act_params_t *act_params)
{
    int i;
    int b_index;
    kappa = act_params->kappa;
    beta = act_params->beta;
    npf = hmc_params->npf;
    mus = hmc_params->mu;
    res_act = hmc_params->res_act;
    res_frc = hmc_params->res_frc;
    eps = res_frc;
    nmax = 1e3;
    
    #if MEASURE_DH
        sprintf(filename_dH, "%s%s", DATA_FOLDER, filename_dH_);
        file_dH = fopen(filename_dH, "w");
        fprintf(file_dH, "time,dH\n");
        fclose(file_dH);
    #endif
    #if MEASURE_ACC
        sprintf(filename_acc, "%s%s", DATA_FOLDER, filename_acc_);
        file_acc = fopen(filename_acc, "w");
        fprintf(file_acc, "time,acc\n");
        fclose(file_acc);
    #endif
    #if MEASURE_PLAQ
        sprintf(filename_plaq, "%s%s", DATA_FOLDER, filename_plaq_);
        file_plaq = fopen(filename_plaq, "w");
        fprintf(file_plaq, "time,plaq\n");
        fclose(file_plaq);
    #endif
    #if MEASURE_QTOP    
        sprintf(filename_qtop, "%s%s", DATA_FOLDER, filename_qtop_);
        file_qtop = fopen(filename_qtop, "w");
        fprintf(file_qtop, "time,qtop\n");
        fclose(file_qtop);
    #endif    
    #if MEASURE_CORR
        /* allocate workspace */
        phi = malloc(2*sizeof(spinor*));
        for (b_index=0; b_index<2; b_index++) {
            phi[b_index] = malloc(V*sizeof(spinor));
        }
        eta = malloc(V*sizeof(spinor));
        /* PP */
        sprintf(filename_PP, "%s%s", DATA_FOLDER, filename_PP_);
        file_PP = fopen(filename_PP, "w");
        fprintf(file_PP, "time,");
        for (i=0; i<V; i++) { 
            fprintf(file_PP, "%d", i);
            if (i != V-1) { fprintf(file_PP, ",");}
        }
        fprintf(file_PP, "\n");
        fclose(file_PP);
        /* A0P */
        sprintf(filename_A0P, "%s%s", DATA_FOLDER, filename_A0P_);
        file_A0P = fopen(filename_A0P, "w");
        fprintf(file_A0P, "time,");
        for (i=0; i<V; i++) { 
            fprintf(file_A0P, "%d", i);
            if (i != V-1) { fprintf(file_A0P, ",");}
        }
        fprintf(file_A0P, "\n");
        fclose(file_A0P);
        /* A1P */
        sprintf(filename_A1P, "%s%s", DATA_FOLDER, filename_A1P_);
        file_A1P = fopen(filename_A1P, "w");
        fprintf(file_A1P, "time,");
        for (i=0; i<V; i++) { 
            fprintf(file_A1P, "%d", i);
            if (i != V-1) { fprintf(file_A1P, ",");}
        }
        fprintf(file_A1P, "\n");
        fclose(file_A1P);
        /* V0V0 */
        sprintf(filename_V0V0, "%s%s", DATA_FOLDER, filename_V0V0_);
        file_V0V0 = fopen(filename_V0V0, "w");
        fprintf(file_V0V0, "time,");
        for (i=0; i<V; i++) { 
            fprintf(file_V0V0, "%d", i);
            if (i != V-1) { fprintf(file_V0V0, ",");}
        }
        fprintf(file_V0V0, "\n");
        fclose(file_V0V0);
        /* V0V1 */
        sprintf(filename_V0V1, "%s%s", DATA_FOLDER, filename_V0V1_);
        file_V0V1 = fopen(filename_V0V1, "w");
        fprintf(file_V0V1, "time,");
        for (i=0; i<V; i++) { 
            fprintf(file_V0V1, "%d", i);
            if (i != V-1) { fprintf(file_V0V1, ",");}
        }
        fprintf(file_V0V1, "\n");
        fclose(file_V0V1);
        /* V1V0 */
        sprintf(filename_V1V0, "%s%s", DATA_FOLDER, filename_V1V0_);
        file_V1V0 = fopen(filename_V1V0, "w");
        fprintf(file_V1V0, "time,");
        for (i=0; i<V; i++) { 
            fprintf(file_V1V0, "%d", i);
            if (i != V-1) { fprintf(file_V1V0, ",");}
        }
        fprintf(file_V1V0, "\n");
        fclose(file_V1V0);
        /* V1V1 */
        sprintf(filename_V1V1, "%s%s", DATA_FOLDER, filename_V1V1_);
        file_V1V1 = fopen(filename_V1V1, "w");
        fprintf(file_V1V1, "time,");
        for (i=0; i<V; i++) { 
            fprintf(file_V1V1, "%d", i);
            if (i != V-1) { fprintf(file_V1V1, ",");}
        }
        fprintf(file_V1V1, "\n");
        fclose(file_V1V1);
    #endif
}

void measure_and_record(spinor** pf, double time, double dH, int acc)
{

    int i;
    /*
    compute_2pt();
    */

    #if MEASURE_DH
        file_dH = fopen(filename_dH, "a");
        fprintf(file_dH, "%4.3e,%4.3e\n", time, dH);
        fclose(file_dH);
    #endif
    #if MEASURE_ACC
        file_acc = fopen(filename_acc, "a");
        fprintf(file_acc, "%4.3e,%i\n", time, acc);
        fclose(file_acc);
    #endif
    #if MEASURE_PLAQ
        plaq = gauge_action(-1)/V;
        file_plaq = fopen(filename_plaq, "a");
        fprintf(file_plaq, "%4.3e,%4.3e\n", time, plaq);
        fclose(file_plaq);
    #endif
    #if MEASURE_QTOP    
        qtop = top_charge();
        file_qtop = fopen(filename_qtop, "a");
        fprintf(file_qtop, "%4.3e,%4.3e\n", time, qtop);
        fclose(file_qtop);
    #endif    
    #if MEASURE_CORR
        compute_correlators();
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
    #endif
}

void run_plot_scripts()
{
    char *c_dH = "";
    char *c_acc = "";
    char *c_plaq = "";
    char *c_qtop = "";
    char *c_corr = "";
    char command[256];
    int out;
    #if MEASURE_DH
        c_dH = "dH";
    #endif
    #if MEASURE_ACC
        c_acc = "acc";
    #endif
    #if MEASURE_PLAQ
        c_plaq = "plaq";
    #endif
    #if MEASURE_QTOP
        c_qtop = "qtop";
    #endif
    #if MEASURE_CORR
        c_corr = "corr";
    #endif
    sprintf(command, "cd ../results/data_analysis/ && python plot.py %i %i %4.3e %4.3e %s %s %s %s %s", T, L, kappa, beta, c_dH, c_acc, c_plaq, c_qtop, c_corr);
    out = system("");
    out = system(command);
    if (out == 0){
        printf("plotting successful\n");
    }
}