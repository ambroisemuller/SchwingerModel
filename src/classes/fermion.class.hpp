/**
 * @file    fermion.class.hpp
 * @brief   Pseudofermion field(s) for the Schwinger model.
 *          Lattice point corresponding to coordinates (t, x) are at index [t*L + x].
 * 
 * @author  Ambroise Muller
 * @date    November 2023
*/

class PseudoFermionField : public Field<Spinor*> {

private:

    bool need_workspace;            // true for physical fields, false for workspace fields
    PseudoFermionField *v1, *v2;    // general workspace fields

    bool need_cg_workspace;
    PseudoFermionField *cg_work, *cg_p, *cg_res, *cg_ap;    // CG workspace

    HoppingField    *hop;           // Hopping field
    GaugeField      *gauge;         // Gauge field

public:

    /**
     * @brief Constructor for PseudoFermionField.
     * @param T_ Temporal extent of lattice.
     * @param L_ Spatial extent of lattice. 
     * @param need_workspace_ true for physical (HMC) fields, false for workspace fields.
    */
    PseudoFermionField(int T_, int L_, bool need_workspace_, bool need_cg_workspace_) : 
    need_workspace(need_workspace_), need_cg_workspace(need_cg_workspace_), Field(T_, L_) {
        for (int i=0; i<V; i++){
            values[i] = new Spinor();
        }
        if (need_workspace_){
            allocate_workspace();
        }
        if (need_cg_workspace_){
            allocate_cg_workspace();
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
            if (need_cg_workspace){
                delete cg_work->values[i];
                delete cg_p->values[i];
                delete cg_res->values[i];
                delete cg_ap->values[i];
            }
        }
        if (need_workspace){
            delete[] v1->values;
            delete[] v2->values;
            delete v1;
            delete v2;
        }
        if (need_cg_workspace){
            delete[] cg_work->values;
            delete[] cg_p->values;
            delete[] cg_res->values;
            delete[] cg_ap->values;
            delete cg_work;
            delete cg_p;
            delete cg_res;
            delete cg_ap;
        }
    }

    /**
     * @brief Allocate memory for auxilliary field instances used in computations.
    */
    void allocate_workspace(){
        v1 = new PseudoFermionField(T, L, false, true);
        v2 = new PseudoFermionField(T, L, false, true);
        for (int i=0; i<V; i++){
            v1->values[i] = new Spinor();
            v2->values[i] = new Spinor();
        }
    }

    void allocate_cg_workspace(){
        cg_work = new PseudoFermionField(T, L, false, false);
        cg_p    = new PseudoFermionField(T, L, false, false);
        cg_res  = new PseudoFermionField(T, L, false, false);
        cg_ap   = new PseudoFermionField(T, L, false, false);
        for (int i=0; i<V; i++){
            cg_work->values[i] = new Spinor();
            cg_p->values[i]    = new Spinor();
            cg_res->values[i]  = new Spinor();
            cg_ap->values[i]   = new Spinor();
        }
    }

    /**
     * @brief Assign hopping field.
     * @param hopping Pointer to hopping field.
    */
    void assign_hopping_field(HoppingField *hopping){
        hop = hopping;
        if (need_workspace){
            // v1->hop = hopping;
            v1->assign_hopping_field(hopping);
            // v2->hop = hopping;
            v2->assign_hopping_field(hopping);
        }
        if (need_cg_workspace){
            cg_work->hop = hopping;
            cg_p->hop    = hopping;
            cg_res->hop  = hopping;
            cg_ap->hop   = hopping;
        }
        Log::print("Successfully assigned hopping field", Log::VERBOSE);
    }

    /**
     * @brief Assign gauge field.
     * @param gauge_ Pointer to gauge field.
    */
    void assign_gauge_field(GaugeField *gauge_){
        gauge = gauge_;
        if (need_workspace){
            // v1->gauge = gauge_;
            v1->assign_gauge_field(gauge_);
            // v2->gauge = gauge_;
            v2->assign_gauge_field(gauge_);
        }
        if (need_cg_workspace){
            cg_work->gauge = gauge_;
            cg_p->gauge    = gauge_;
            cg_res->gauge  = gauge_;
            cg_ap->gauge   = gauge_;
        }
        Log::print("Successfully assigned gauge field", Log::VERBOSE);
    }

    void set_to_zero(){
        for (int i=0; i<V; i++){
            values[i]->s[0] = 0;
            values[i]->s[1] = 0;
        }
    }

    void set_to(PseudoFermionField *source){
        for (int i=0; i<V; i++){
            values[i]->s[0] = source->values[i]->s[0];
            values[i]->s[1] = source->values[i]->s[1];
        }
    }

    void apply_gamma5(){
        for (int i=0; i<V; i++){
            values[i]->s[1] = -values[i]->s[1];
        }
    }

    void add_z_times_pf(complex<double> z, PseudoFermionField *pf){
        for (int i=0; i<V; i++){
            values[i]->s[0] += z * pf->values[i]->s[0];
            values[i]->s[1] += z * pf->values[i]->s[1];
        }
    }

    void set_to_pf1_plus_z_times_pf2(PseudoFermionField *pf1, complex<double> z, PseudoFermionField *pf2){
        for (int i=0; i<V; i++){
            values[i]->s[0] = pf1->values[i]->s[0] + z * pf2->values[i]->s[0];
            values[i]->s[1] = pf1->values[i]->s[1] + z * pf2->values[i]->s[1];
        }
    }

    /**
     * @brief Computes the scalar product with another spinor field.
     * @param other Another spinor/pseudofermion field.
     * @return Scalar product = sum_{sites i} |conj(this_i) * other_i|^2.
    */
    complex<double> scalar_prod_with(PseudoFermionField *other){
        complex<double> res = 0.0;
        for (int i=0; i<V; i++){
            res += conj(values[i]->s[0])*other->values[i]->s[0];
            res += conj(values[i]->s[1])*other->values[i]->s[1];
        }
        return res;
    }


    double compute_fermion_action1(double kappa, double mu, double res){
        v1->set_to_CG_solution(this, kappa, mu, res, N_MAX);
        return scalar_prod_with(v1).real();
    }

    double compute_fermion_action2(double kappa, double mu1, double mu2, double res){
        return (mu2*mu2 - mu1*mu1)*compute_fermion_action1(kappa, mu1, res);
    }

    double setpf1(double kappa, double mu){
        Random::gauss_spinor_field(V, v1->values);
        set_to_dirac_of(v1, kappa, mu);
        apply_gamma5();
        return v1->scalar_prod_with(v1).real();
    }

    double setpf2(double kappa, double mu1, double mu2, double res){
        Random::gauss_spinor_field(V, values);
        v1->set_to_CG_solution(this, kappa, mu2, res, N_MAX);
        v2->set_to(v1);
        v1->set_to_dirac_of(v2, kappa, -mu2);
        v1->apply_gamma5();
        add_z_times_pf(COMPLEX_I*(mu1-mu2), v1);
        return (mu2*mu2 - mu1*mu1) * v1->scalar_prod_with(v1).real();
    }

    void compute_force1_to(ForceField *force, double kappa, double mu, double a, double res){
        v1->set_to_CG_solution(this, kappa, mu, res, N_MAX);
        v2->set_to_dirac_of(v1, kappa, mu);
        v2->apply_gamma5();
        for (int i=0; i<V; i++){
            for (int nu=0; nu<D; nu++){
                complex<double> c = -a * kappa * exp(- COMPLEX_I * gauge->values[i][nu]);
                Spinor tmp;
                mul_1_plus_gamma(nu, v2->values[hop->values[i][nu]], &tmp);
                force->values[i][nu] = 2 * ((conj(v1->values[i]->s[0])*c*tmp.s[0]).imag() - (conj(v1->values[i]->s[1])*c*tmp.s[1]).imag());
                mul_1_plus_gamma(nu, v1->values[hop->values[i][nu]], &tmp);
                force->values[i][nu] = 2 * ((conj(v2->values[i]->s[0])*c*tmp.s[0]).imag() - (conj(v2->values[i]->s[1])*c*tmp.s[1]).imag());
            }
        }

    }

    void compute_force2_to(ForceField *force, double kappa, double mu1, double mu2, double a, double res){
        compute_force1_to(force, kappa, mu1, a*(mu2*mu2 - mu1*mu1), res);
    }

    /**
     * @brief Conjugate gradient algorithm to find solution to (D^dag*D + mu^2) (*this) = (*other).
     *        Solution is written to values of current instance.
     * @param other Pointer to PseudoFermionField instance on RHS of equation to solve.
     * @param kappa Mass parameter.
     * @param mu    Twisted mass.
     * @param eps   Solver tolerance.
     * @param nmax  Maximum number of solver iterations.
    */
    void set_to_CG_solution(PseudoFermionField *other, double kappa, double mu, double eps, int nmax){
        set_to_zero();
        cg_res->set_to(other);
        cg_p->set_to(cg_res);
        double rsold = cg_res->scalar_prod_with(cg_res).real();
        for (int i=0; i<nmax; i++){
            cg_work->set_to_dirac_of(cg_p, kappa, mu);
            cg_work->apply_gamma5();
            cg_ap->set_to_dirac_of(cg_work, kappa, -mu);
            cg_ap->apply_gamma5();
            complex<double> c = cg_p->scalar_prod_with(cg_ap);
            double alpha = rsold/c.real();
            add_z_times_pf(alpha, cg_p);
            cg_res->add_z_times_pf(-alpha, cg_ap);
            double rsnew = cg_res->scalar_prod_with(cg_res).real();
            if (sqrt(rsnew) < eps){
                break;
            }
            cg_p->set_to_pf1_plus_z_times_pf2(cg_res, rsnew/rsold, cg_p);
            rsold = rsnew;
            if (i == nmax-1){
                Log::print("Error: CG did not converge", Log::ERROR);
                break; // @todo remove
                exit(1);
            }
        }
    }

    /**
     * @brief Set values of current field to D*other, where D is the Dirac operator.
     * @param other Source PseudoFermionField instance.
     * @param kappa Mass parameter.
     * @param mu    Twisted mass.
    */
    void set_to_dirac_of(PseudoFermionField *other, double kappa, double mu){
        for (int i=0; i<V; i++){
            values[i]->s[0] = complex<double>(1, mu) * other->values[i]->s[0];
            values[i]->s[1] = complex<double>(1, -mu) * other->values[i]->s[1];
            complex<double> c;
            Spinor tmp;
            for (int nu=0; nu<D; nu++){
                if (nu < D){    // forward hopping
                    c = -kappa * exp(-COMPLEX_I * gauge->values[i][nu]);
                } else {        // backward hopping
                    c = -kappa * exp(COMPLEX_I * gauge->values[hop->values[i][nu]][nu-D]);
                }
                mul_1_plus_gamma(nu, other->values[hop->values[i][nu]], &tmp);
                values[i]->s[0] += c * tmp.s[0];
                values[i]->s[1] += c * tmp.s[1];
            }
        }
    }

};
