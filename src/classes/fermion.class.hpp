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

    bool need_workspace;            // true for physical fields, false for workspace fields.
    PseudoFermionField *v1, *v2;    // workspace fields

    HoppingField    *hop;           // Hopping field

    /**
     * @brief Constructor for PseudoFermionField.
     * @param T_ Temporal extent of lattice.
     * @param L_ Spatial extent of lattice. 
     * @param need_workspace_ true for physical (HMC) fields, false for workspace fields.
    */
    PseudoFermionField(int T_, int L_, bool need_workspace_) : 
    need_workspace(need_workspace_), Field(T_, L_) {
        for (int i=0; i<V; i++){
            values[i] = new Spinor();
        }
        if (need_workspace_){
            allocate_workspace();
        }
    }

    /**
     * @brief Destructor for PseudoFermionField.
    */
    ~PseudoFermionField() {
        for (int i=0; i<V; i++){
            delete values[i];
            if (need_workspace){
                delete v1->values[i];
                delete v2->values[i];
            }
        }
        if (need_workspace){
            delete[] v1->values;
            delete v1;
            delete[] v2->values;
            delete v2;
        }
    }

    /**
     * @brief Allocate memory for auxilliary field instances used in computations.
    */
    void allocate_workspace(){
        v1 = new PseudoFermionField(T, L, false);
        v2 = new PseudoFermionField(T, L, false);
        for (int i=0; i<V; i++){
            v1->values[i] = new Spinor();
            v2->values[i] = new Spinor();
        }
    }

    /**
     * @brief Assign hopping field.
     * @param hopping Pointer to hopping field.
    */
    void assign_hopping_field(HoppingField *hopping){
        hop = hopping;
        Log::print("Successfully assigned hopping field", Log::VERBOSE);
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

    /**
     * @brief Conjugate gradient algorithm to find solution to (D^dag*D + mu^2) (*this) = (*other).
     *        Solution is written to values of current instance.
     * @param other Pointer to PseudoFermionField instance on RHS of equation to solve.
     * @param kappa 
    */
    void set_to_CG_solution(PseudoFermionField *other, double kappa, double mu){

    }

};