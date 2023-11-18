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

    const int   D = 2;      // Dimension (fixed)
    const int   T;          // Number of sites in temporal direction
    const int   L;          // Number of sites in spatial direction
    const int   V;          // Lattice volume
    double      a;          // Lattice spacing

    GaugeField  *gauge;      // Gauge field
    int         **hop;      // Hopping field

    Lattice(int T_, int L_, double a_) : T(T_), L(L_), V(T_*L_), a(a_) {

        Log::print("Building lattice of size "+to_string(T_)+" by "+to_string(L_)+" with spacing "+to_string(a), Log::VERBOSE);

        gauge = new GaugeField(T, L);
        gauge->initialize_cold();

        hop = new int*[V];

        for (int i=0; i<V; i++){
            hop[i] = new int[2*D];
        }

        Log::print("Successfully initialized gauge and hopping fields", Log::VERBOSE);

    }

};