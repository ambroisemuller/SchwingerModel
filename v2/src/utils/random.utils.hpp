/**
 * @file    random.utils.hpp
 * @brief   Random number generator.
 *          Wrapper of ranlxs/ranlxd algorithm developed by Martin Luescher.
 * 
 * @author  Ambroise Muller
 * @date    November 2023
*/


#define BASE 0x1000000
#define MASK 0xffffff

#define STEP(pi, pj) \
        d = (*pj).c1.c1 - (*pi).c1.c1 - carry.c1; \
        (*pi).c2.c1 += (d < 0); \
        d += BASE; \
        (*pi).c1.c1 = d&MASK; \
        d = (*pj).c1.c2 - (*pi).c1.c2 - carry.c2; \
        (*pi).c2.c2 += (d < 0); \
        d += BASE; \
        (*pi).c1.c2 = d&MASK; \
        d = (*pj).c1.c3 - (*pi).c1.c3 - carry.c3; \
        (*pi).c2.c3 += (d < 0); \
        d += BASE; \
        (*pi).c1.c3 = d&MASK; \
        d = (*pj).c1.c4 - (*pi).c1.c4 - carry.c4; \
        (*pi).c2.c4 += (d < 0); \
        d += BASE; \
        (*pi).c1.c4 = d&MASK; \
        d = (*pj).c2.c1 - (*pi).c2.c1; \
        carry.c1 = (d < 0); \
        d += BASE; \
        (*pi).c2.c1 = d&MASK; \
        d = (*pj).c2.c2 - (*pi).c2.c2; \
        carry.c2 = (d < 0); \
        d += BASE; \
        (*pi).c2.c2 = d&MASK; \
        d = (*pj).c2.c3 - (*pi).c2.c3; \
        carry.c3 = (d < 0); \
        d += BASE; \
        (*pi).c2.c3 = d&MASK; \
        d = (*pj).c2.c4 - (*pi).c2.c4; \
        carry.c4 = (d < 0); \
        d += BASE; \
        (*pi).c2.c4 = d&MASK


typedef struct
{
   int c1, c2, c3, c4;
} vec_t;

typedef struct
{
   vec_t c1, c2;
} dble_vec_t;


namespace Random {

    static int init = 0, pr, prm, ir, jr, is, is_old, next[96];
    static double one_bit;
    static vec_t carry;

    static union{
        dble_vec_t vec[12];
        int num[96];
    } x;

    static void error(int no){
        switch(no)
        {
            case 0:
                Log::print("Error in rlxd_init\nArithmetic on this machine is not suitable for ranlxd", Log::ERROR);
                break;
            case 1:
                Log::print("Error in subroutine rlxd_init\nBad choice of luxury level (should be 1 or 2)", Log::ERROR);
                break;
            case 2:
                Log::print("Error in subroutine rlxd_init\nBad choice of seed (should be between 1 and 2^31-1)", Log::ERROR);
                break;
            case 3:
                Log::print("Error in rlxd_get\nUndefined state (ranlxd is not initialized)", Log::ERROR);
                break;
            case 4:
                Log::print("Error in rlxd_reset\nArithmetic on this machine is not suitable for ranlxd", Log::ERROR);
                break;
            case 5:
                Log::print("Error in rlxd_reset\nUnexpected input data", Log::ERROR);
                break;
        }      
        Log::print("Program aborted", Log::ERROR);
        exit(0);
    };


    static void update(){
        int k, kmax, d;
        dble_vec_t *pmin, *pmax, *pi, *pj;
        kmax = pr;
        pmin = &x.vec[0];
        pmax = pmin + 12;
        pi = &x.vec[ir];
        pj = &x.vec[jr];
        for (k=0; k<kmax; k++){
            STEP(pi, pj);
            pi += 1;
            pj += 1;
            if (pi == pmax){
                pi = pmin;
            }
            if (pj == pmax){
                pj = pmin;
            } 
        }
        ir += prm;
        jr += prm;
        if (ir >= 12){
            ir -= 12;
        }
        if (jr >= 12){
            jr -= 12;
        }
        is = 8 * ir;
        is_old = is;
    }


    static void define_constants(){
        int k;
        one_bit = ldexp(1.0,-24);
        for (k=0; k<96; k++){
            next[k] = (k+1) % 96;
            if ((k%4) == 3) {
                next[k] = (k+5) % 96;
            }
        }   
    }


    void rlxd_init(int level, int seed){
        int i, k, l;
        int ibit, jbit, xbit[31];
        int ix, iy;
        if ((INT_MAX<2147483647) || (FLT_RADIX!=2) || (FLT_MANT_DIG<24) || (DBL_MANT_DIG<48)) {
            error(0);
        }
        define_constants();
        if (level == 1){
            pr = 202;
        }
        else if (level == 2){
            pr = 397;
        }
        else {
            error(1);
        }
        i = seed;
        for (k=0; k<31; k++){
            xbit[k] = i % 2;
            i /= 2;
        }
        if ((seed<=0) || (i!=0)){
            error(2);
        }
        ibit = 0;
        jbit = 18;
        for (i=0; i<4; i++){
            for (k=0; k<24; k++){
                ix = 0;
                for (l=0; l<24; l++){
                    iy = xbit[ibit];
                    ix = 2 * ix + iy;
                    xbit[ibit] = (xbit[ibit] + xbit[jbit]) % 2;
                    ibit = (ibit+1) % 31;
                    jbit = (jbit+1) % 31;
                }
                if ((k%4) != i){
                    ix = 16777215 - ix;
                }
                x.num[4*k+i] = ix;
            }
        }
        carry.c1 = 0;
        carry.c2 = 0;
        carry.c3 = 0;
        carry.c4 = 0;
        ir = 0;
        jr = 7;
        is = 91;
        is_old = 0;
        prm = pr % 12;
        init = 1;
    }


    void ranlxd(double r[],int n){
        int k;
        if (init == 0){
            rlxd_init(1, 1);
        }
        for (k=0; k<n; k++){
            is = next[is];
            if (is == is_old){
                update();
            }
            r[k] = one_bit * ((double)(x.num[is+4]) + one_bit*(double)(x.num[is]));      
        }
    }


    int rlxd_size(){
        return(105);
    }


    void rlxd_get(int state[]){
        int k;
        if (init == 0){
            error(3);
        }
        state[0] = rlxd_size();
        for (k=0; k<96; k++)
            state[k+1] = x.num[k];
        state[97] = carry.c1;
        state[98] = carry.c2;
        state[99] = carry.c3;
        state[100] = carry.c4;
        state[101] = pr;
        state[102] = ir;
        state[103] = jr;
        state[104] = is;
    }


    void rlxd_reset(int state[]){
        int k;
        if ((INT_MAX<2147483647) || (FLT_RADIX!=2) || (FLT_MANT_DIG<24) || (DBL_MANT_DIG<48)){
            error(4);
        }
        define_constants();
        if (state[0] != rlxd_size()){
            error(5);
        }
        for (k=0;k<96;k++){
            if ((state[k+1]<0) || (state[k+1]>=167777216)){
                error(5);
            }
            x.num[k] = state[k+1];
        }
        if (((state[97]!=0)&&(state[97]!=1)) || ((state[98]!=0)&&(state[98]!=1)) || ((state[99]!=0)&&(state[99]!=1)) || ((state[100]!=0)&&(state[100]!=1))){
            error(5);
        }
        carry.c1 = state[97];
        carry.c2 = state[98];
        carry.c3 = state[99];
        carry.c4 = state[100];
        pr = state[101];
        ir = state[102];
        jr = state[103];
        is = state[104];
        is_old = 8 * ir;
        prm = pr % 12;
        init = 1;
        if (((pr!=202)&&(pr!=397)) || (ir<0) || (ir>11) || (jr<0) || (jr>11) || (jr!=((ir+7)%12)) || (is<0)||(is>91)){
            error(5);
        }
    }

    double gauss_normalization = 8.0*atan(1.0);

    /**
     * @brief Initialize array to normally distributed random real values.
     * @param n Length of doubles array.
     * @param rand Array of doubles to initialize.
    */
    void gauss_rand(int n, double *rand){
        int i = 0;
        double tmp[2], r, phi;
        while (i < n){
            ranlxd(tmp, 2); // two random numbers, flat distribution in [0,1)
            phi = tmp[0] * gauss_normalization; // compute polar coordinates: angle and radius
            r = sqrt(-log(1.0-tmp[1])); // map second number [0,1) -> (0,1]
            rand[i] = r * cos(phi); 
            i++;
            if (i < n) {    // compute second only if requested
                rand[i] = r * sin(phi);
            }; 
            i++;
        }
    }

    /**
     * @brief Initialize array to normally distributed random real values.
     * @param n Length of doubles array.
     * @param rand Array of pointers to doubles to initialize.
    */
    void gauss_rand(int n, double **rand){
        int i = 0;
        double tmp[2], r, phi;
        while (i < n){
            ranlxd(tmp, 2); // two random numbers, flat distribution in [0,1)
            phi = tmp[0] * gauss_normalization; // compute polar coordinates: angle and radius
            r = sqrt(-log(1.0-tmp[1])); // map second number [0,1) -> (0,1]
            *(rand[i]) = r * cos(phi); 
            i++;
            if (i < n) {    // compute second only if requested
                *(rand[i]) = r * sin(phi);
            }; 
            i++;
        }
    }

    /**
     * @brief Initialize array to normally distributed random spinor values.
     * @param field_volume Length of spinors array.
     * @param field_values Array of pointers to spinor values to initialize.
    */
    void gauss_spinor_field(int field_volume, Spinor **field_values){
        double rv[4];
        for (int i=0; i<field_volume; i++){
            gauss_rand(4, rv);
            field_values[i]->s[0] = complex<double>(rv[0], rv[1]);
            field_values[i]->s[1] = complex<double>(rv[2], rv[3]);
        }
    }

};
