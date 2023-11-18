/**
 * @file    lattice.class.hpp
 * @brief   Definition of principal Lattice class.
 *          Lattice point corresponding to coordinates (t, x) are at index [t*L + x].
 * 
 * @author  Ambroise Muller
 * @date    November 2023
*/

#include "hopping.class.hpp"
#include "gauge.class.hpp"

class Lattice {

public:

    const int   T;          // Number of sites in temporal direction
    const int   L;          // Number of sites in spatial direction
    const int   V;          // Lattice volume
    double      a;          // Lattice spacing

    GaugeField      *gauge;     // Gauge field
    HoppingField    *hop;       // Hopping field

private:

    double          time;           // Simulation time
    GaugeField      *gauge_old;     // Gauge field from previous step

    #if MEASURE_DH
        double dH;
        stringstream* csv_dH = new std::stringstream;
    #endif
    #if MEASURE_ACC
        int acc;
        stringstream* csv_acc = new std::stringstream;
    #endif
    #if MEASURE_PLAQUETTE
        double plaq;
        stringstream* csv_plaq = new std::stringstream;
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

        hop = new HoppingField(T_, L_);
        hop->initialize();

        Log::print("Successfully initialized gauge and hopping fields", Log::VERBOSE);

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
        delete hop;
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
        // @todo integrate();
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


    double start_hmc(){
        double H = 0.0;
        // @todo H += mom_heat();
        // @todo H += fermion_heat();
        H += gauge->compute_gauge_action(BETA, hop);
        return H;
    }

    double hamiltonian(){
        return -79;
    }

    void record_observables() {
        #if MEASURE_DH
            *csv_dH << time << "," << dH << endl;
        #endif
        #if MEASURE_ACC
            *csv_acc << time << "," << acc << endl;
        #endif
        #if MEASURE_PLAQUETTE
            plaq = gauge->compute_gauge_action(-1, hop)/V;
            *csv_plaq << time << "," << plaq << endl;
        #endif
        #if MEASURE_QTOP
            qtop = gauge->compute_topological_charge(hop);
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