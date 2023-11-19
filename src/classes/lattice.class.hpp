/**
 * @file    lattice.class.hpp
 * @brief   Definition of principal Lattice class.
 *          Lattice point corresponding to coordinates (t, x) are at index [t*L + x].
 * 
 * @author  Ambroise Muller
 * @date    November 2023
*/

#include "field.class.hpp"
#include "force.class.hpp"
#include "momentum.class.hpp"
#include "gauge.class.hpp"
#include "fermion.class.hpp"

class Lattice {

public:

    const int   T;                      // Number of sites in temporal direction
    const int   L;                      // Number of sites in spatial direction
    const int   V;                      // Lattice volume
    double      a;                      // Lattice spacing

    GaugeField          *gauge;         // Gauge field
    PseudoFermionField  **pfermion;     // Pseudofermion fields

    double mu_list[N_PF] = MU_LIST;     // Hasenbusch masses

private:

    double  time;                       // Simulation time
    double  eps = TRAJ_LENGTH/N_STEP;   // Molecular dynamics step size  

    GaugeField      *gauge_old;         // Gauge field from previous step
    MomentumField   *momentum;          // Momentum field
    ForceField      *force;             // Force field
    ForceField      *force_tmp;         // Auxilliary force field (workspace)

    double dH;                          // Change in energy over trajectory
    double acc;                         // Metropolis acceptance flag

    #if MEASURE_DH
        stringstream* csv_dH = new stringstream;
    #endif
    #if MEASURE_ACC
        stringstream* csv_acc = new stringstream;
    #endif
    #if MEASURE_PLAQUETTE
        double plaq;
        stringstream* csv_plaq = new stringstream;
    #endif
    #if MEASURE_QTOP
        double qtop;
        stringstream* csv_qtop = new std::stringstream;
    #endif
    #if MEASURE_PSCC
        double *pscc;
        stringstream* csv_pscc = new std::stringstream;
    #endif

public:

    /**
     * @brief Constructor for Lattice class.
     * @param T_ Temporal extent of lattice.
     * @param L_ Spatial extent of lattice. 
     * @param a_ Lattice spacing.
    */
    Lattice(int T_, int L_, double a_) : T(T_), L(L_), V(T_*L_), a(a_) {

        Log::print("Building lattice of size "+to_string(T_)+" by "+to_string(L_)+" with spacing "+to_string(a), Log::VERBOSE);

        gauge = new GaugeField(T_, L_);
        gauge->initialize_cold();
        gauge_old = new GaugeField(T_, L_);
        gauge_old->initialize_cold();
        Log::print("Successfully initialized gauge fields", Log::VERBOSE);

        pfermion = new PseudoFermionField*[N_PF];
        for (int i=0; i<N_PF; i++){
            pfermion[i] = new PseudoFermionField(T_, L_);
        }
        Log::print("Successfully initialized pseudofermion fields", Log::VERBOSE);

        momentum = new MomentumField(T, L);
        Log::print("Successfully initialized momentum field", Log::VERBOSE);

        force = new ForceField(T, L);
        force_tmp = new ForceField(T, L);
        Log::print("Successfully initialized force fields", Log::VERBOSE);

        #if MEASURE_PSCC
            pscc = new double[V];
        #endif

    }

    /**
     * @brief Destructor for Lattice class.
    */
    ~Lattice(){
        delete gauge;
        delete gauge_old;
        for (int i=0; i<N_PF; i++){
            delete pfermion[i];
        }
        delete[] pfermion;
        delete momentum;
        delete force;
        delete force_tmp;
    }

    /**
     * @brief Run Hamiltonian Monte Carlo (HMC) simulation.
    */
    void HMC() {
        
        time = 0;
        auto start = chrono::high_resolution_clock::now();

        Log::progressBar(0);
        for (int itraj=0; itraj<N_TRAJ; itraj++) {
            single_trajectory();
            time += TRAJ_LENGTH;
            record_observables();
            Log::progressBar(itraj/double(N_TRAJ));
        }
        Log::progressBar(1);
        Log::newLine();

        auto end = chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end-start);
        Log::print("HMC simulation completed in "+to_string(duration.count())+" seconds.", Log::GOOD_NEWS);

    }

    /**
     * @brief Simulate a single HMC trajectory.
    */
    void single_trajectory() {
        gauge_old->assign_link(gauge);

        double H1 = start_hmc();
        integrate();
        double H2 = hamiltonian();

        dH = H2 - H1;
        acc = 0;
        
        if (dH < 0) {
            acc = 1;
        } else {
            double r;
            Random::ranlxd(&r, 1);
            if (exp(-dH) > r) { 
                acc=1; 
            } else { 
                gauge->assign_link(gauge_old); 
            }
        }

    }

    /**
     * @brief Initialize all HMC and compute initial Hamiltonian.
     * @return Hamiltonian value before MD trajectory.
    */
    double start_hmc(){
        double H = 0.0;
        H += momentum->initialize_momentum();
        H += fermion_heat();
        H += gauge->compute_gauge_action(BETA);
        return H;
    }

    /**
     * @brief Calls the heatbath for the various pseudofermion fields, accumulating the actions.
     * @return Pseudofermion contribution to the action.
    */
    double fermion_heat(){
        double action = 0.0;
        if (N_PF > 0){
            for (int i=0; i<N_PF-1; i++){
                action += pfermion[i]->setpf2(KAPPA, mu_list[i], mu_list[i+1], RES_ACT);
            }
            action += pfermion[N_PF-1]->setpf1(KAPPA, mu_list[N_PF-1]);
        }   
        return action;
    }

    /**
     * @brief Compute Hamiltonian of HMC evolution, adding various action terms.
     * @return Hamiltonian after MD trajectory.
    */
    double hamiltonian(){
        double action = momentum->compute_momentum_action();
        if (N_PF > 0){   
            for (int i=0; i<N_PF-1; i++){
                action += pfermion[i]->compute_fermion_action2(KAPPA, mu_list[i], mu_list[i+1], RES_ACT);
            }
            action += pfermion[N_PF-1]->compute_fermion_action1(KAPPA, mu_list[N_PF-1], RES_ACT);
        }
        action += gauge->compute_gauge_action(BETA);
        return action;
    }

    /**
     * @brief Compute forces and update momenta.
     * @param eps_ Step size.
    */
    void move_momentum(double eps_){
        force->set_zero();
        if (N_PF > 0){
            for (int i=0; i<N_PF-1; i++){
                pfermion[i]->compute_force2_to(force_tmp, KAPPA, mu_list[i], mu_list[i+1], 1.0, RES_FRC);
                *force += *force_tmp;
            }
            pfermion[N_PF-1]->compute_force1_to(force_tmp, KAPPA, mu_list[N_PF-1], 1.0, RES_FRC);
            *force += *force_tmp;
        }
        gauge->compute_force_to(force_tmp, BETA, 1.0);
        *force += *force_tmp;
        for (int i=0; i<V; i++){
            for (int nu=0; nu<D; nu++){
                momentum->values[i][nu] -= eps_*force->values[i][nu];
            }
        }
    }

    /**
     * @brief Integration scheme for MD evolution of momentum and gauge fields.
    */
    void integrate(){
        #if INT == LFRG
            for (int i=0; i<N_STEP; i++){
                move_momentum(0.5*eps);
                gauge->move_gauge(eps, momentum);
                move_momentum(0.5*eps);
            }
        #elif INT == OMF2
            for (int i=0; i<N_STEP; i++){
                move_momentum(LAMBDA*eps);
                gauge->move_gauge(0.5*eps, momentum);
                move_momentum((1.0-2.0*LAMBDA)*eps);
                gauge->move_gauge(0.5*eps, momentum);
                move_momentum(LAMBDA*eps);
            }
        #elif INT == OMF4
            double  r1=0.08398315262876693,
                    r2=0.2539785108410595,
                    r3=0.6822365335719091,
                    r4=-0.03230286765269967;
            for (int i=0; i<N_STEP; i++){
                move_momentum(r1*eps);
                gauge->move_gauge(r2*eps, momentum);
                move_mom(r3*eps);
                gauge->move_gauge(r4*eps, momentum);
                move_mom((0.5-r1-r3)*eps);
                gauge->move_gauge((1-2*(r2+r4))*eps, momentum);
                move_mom((0.5-r1-r3)*eps);
                gauge->move_gauge(r4*eps, momentum);
                move_mom(r3*eps);
                gauge->move_gauge(r2*eps, momentum);
                move_mom(r1*eps);
            }
        #else
            Log::print("Unknown integrator (known: LFRG, OML2, OML4)", Log::ERROR);
            exit(1);
        #endif
    }

    /**
     * @brief Measure desired observables and record values.
    */
    void record_observables() {
        #if MEASURE_DH
            *csv_dH << time << "," << dH << endl;
        #endif
        #if MEASURE_ACC
            *csv_acc << time << "," << acc << endl;
        #endif
        #if MEASURE_PLAQUETTE
            plaq = gauge->compute_gauge_action(-1)/V;
            *csv_plaq << time << "," << plaq << endl;
        #endif
        #if MEASURE_QTOP
            qtop = gauge->compute_topological_charge();
            *csv_qtop << time << "," << qtop << endl;
        #endif
        #if MEASURE_PSCC
            for (int i=0; i<V; i++) {
                pscc[i] = 0; // @todo implement PSCC measurement
            }
            *csv_pscc << time;
            for (int i=0; i<V; i++) {
                *csv_pscc << ", " << pscc[i];
            }
            *csv_pscc << endl;
        #endif
    }

    /**
     * @brief Save measured observables to results folder.
     * @param path Path to results directory.
    */
    void save_results(const string& path) {
        ofstream output;
        #if MEASURE_DH
            output.open(path+"/dH.csv");
            output << csv_dH->str();
            output.close();
        #endif
        #if MEASURE_ACC
            output.open(path+"/acc.csv");
            output << csv_acc->str();
            output.close();
        #endif
        #if MEASURE_PLAQUETTE
            output.open(path+"/plaq.csv");
            output << csv_plaq->str();
            output.close();
        #endif
        #if MEASURE_QTOP
            output.open(path+"/qtop.csv");
            output << csv_qtop->str();
            output.close();
        #endif
        #if MEASURE_PSCC
            output.open(path+"/pscc.csv");
            output << csv_pscc->str();
            output.close();
        #endif
    }

};