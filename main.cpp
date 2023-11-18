/**
 * @file    main.hpp
 * @brief   Main project file.
 * 
 * @author  Ambroise Muller
 * @date    November 2023
*/


#include "src/imports.hpp"


int main(int argc, char* argv[]) {

    if (argc == 1){
    
        Random::rlxd_init(1, SEED);
        Lattice* lattice = new Lattice(Nt, Nx, spacing);

        lattice->gauge->initialize_from("../input/config.csv");
        lattice->HMC();
        lattice->gauge->save_config("../results/gauge_config.csv");
        cout << (*(lattice->gauge))[2][0] << endl;

    } else {
        Log::print("Number of arguments not supported!", Log::ERROR);
        exit(1);
    }

    return 0;
    
}