#define     RES_ACT         1.e-10                          // precision solver for action computation: accept/reject step
#define     RES_FRC         1.e-8                           // precision solver for force computation
#define     N_MAX           1e3                             // Maximum number of steps for CG solver

#define     SEED            116                             // seed for random number generator

#define     OUTPUT_CONFIG   "../results/gauge_config.csv"   // path to store output gauge field configuration

#define     MULTI_TD        false                           // use multithreading on CPUs
#define     GPU             false                           // use CUDA 
