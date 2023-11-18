/**
 * @file    lattice.class.hpp
 * @brief   Define Lattice class.
 * 
 * @author  Ambroise Muller
 * @date    November 2023
*/

#include "gauge.class.hpp"
#include "hopping.class.hpp"

class Lattice {

public:

    const int   T;          // Number of sites in temporal direction
    const int   L;          // Number of sites in spatial direction
    const int   V;          // Lattice volume
    double      a;          // Lattice spacing

    GaugeField      *gauge;     // Gauge field
    HoppingField    *hop;       // Hopping field

    Lattice(int T_, int L_, double a_) : T(T_), L(L_), V(T_*L_), a(a_) {

        Log::print("Building lattice of size "+to_string(T_)+" by "+to_string(L_)+" with spacing "+to_string(a), Log::VERBOSE);

        gauge = new GaugeField(T_, L_);
        gauge->initialize_cold();

        hop = new HoppingField(T_, L_);
        hop->initialize();

        Log::print("Successfully initialized gauge and hopping fields", Log::VERBOSE);

    }

    void HMC() {
        
        int max_steps = 100;
        auto start = chrono::high_resolution_clock::now();

        for (int step=0; step<max_steps; step++) {
            Log::progressBar(step/double(max_steps));
        }
        Log::progressBar(1);
        Log::newLine();

        auto end = chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(end-start);
        Log::print("HMC simulation completed in "+to_string(duration.count())+" seconds.", Log::GOOD_NEWS);
        
    }

};