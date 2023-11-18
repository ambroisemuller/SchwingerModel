/**
 * @file    main.hpp
 * @brief   Main project file.
 *          Schwinger model simulations with Wilson fermions. 
 * 
 * @author  Ambroise Muller
 * @date    November 2023
*/

#include "src/imports.hpp"


int main(int argc, char* argv[]) {

    if (argc == 1){
    
        Random::rlxd_init(1, SEED);
        Lattice* lattice = new Lattice(Nt, Nx, spacing);
        // lattice->gauge->initialize_from("../input/config.csv");
        lattice->HMC();
        lattice->gauge->save_config("../results/gauge_config.csv");
        lattice->save_results(RESULTS_PATH);
        Log::print(string("Results saved to folder ")+(RESULTS_PATH), Log::GOOD_NEWS);

    } else {
        Log::print("Number of arguments not supported!", Log::ERROR);
        exit(1);
    }

    return 0;
    
}