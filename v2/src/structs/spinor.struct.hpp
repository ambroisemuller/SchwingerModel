
#define     COMPLEX_I   complex<double>(0, 1)

typedef struct {
    complex<double> s[2];
} Spinor;


void mul_1_plus_gamma(int mu, Spinor* source, Spinor* target){
    if (mu == 0){
        target->s[0] = source->s[0] - source->s[1];
        target->s[1] = - target->s[0];
    }
    else if (mu == 1){
        target->s[0] = source->s[0] + COMPLEX_I * source->s[1];
        target->s[1] = - COMPLEX_I * target->s[0];
    }
    else if (mu == 2){
        target->s[0] = source->s[0] + source->s[1];
        target->s[1] = target->s[0];
    } 
    else if (mu == 3){
        target->s[0] = source->s[0] - COMPLEX_I * source->s[1];
        target->s[1] = COMPLEX_I * target->s[0];
    }
}