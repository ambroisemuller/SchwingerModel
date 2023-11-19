/**
 * @file    fermion.class.hpp
 * @brief   Pseudofermion field(s) for the Schwinger model.
 *          Lattice point corresponding to coordinates (t, x) are at index [t*L + x].
 * 
 * @author  Ambroise Muller
 * @date    November 2023
*/


#include "../structs/spinor.struct.hpp"


class PseudoFermionField : public Field<Spinor*> {

public:

    PseudoFermionField(int T_, int L_) : Field(T_, L_) {
        for (int i=0; i<V; i++){
            values[i] = new Spinor();
        }
    }

    ~PseudoFermionField() {
        for (int i=0; i<V; i++){
            delete values[i];
        }
    }

    double compute_fermion_action1(double kappa, double mu, double res){
        return 0;
    }


    double compute_fermion_action2(double kappa, double mu1, double mu2, double res){
        return 0;
    }

    double setpf1(double kappa, double mu){
        return 0;
    }

    double setpf2(double kappa, double mu1, double mu2, double res){
        return 0;
    }

    void compute_force1_to(ForceField *force, double kappa, double mu, double a, double res){

    }

    void compute_force2_to(ForceField *force, double kappa, double mu1, double mu2, double a, double res){

    }

};