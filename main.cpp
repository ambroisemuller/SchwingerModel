/**
 * @file    main.cpp
 * @brief   Main project file.
 *          Schwinger model simulations with Wilson fermions. 
 * 
░██████╗░█████╗░██╗░░██╗░██╗░░░░░░░██╗██╗███╗░░██╗░██████╗░███████╗██████╗░  ███╗░░░███╗░█████╗░██████╗░███████╗██╗░░░░░
██╔════╝██╔══██╗██║░░██║░██║░░██╗░░██║██║████╗░██║██╔════╝░██╔════╝██╔══██╗  ████╗░████║██╔══██╗██╔══██╗██╔════╝██║░░░░░
╚█████╗░██║░░╚═╝███████║░╚██╗████╗██╔╝██║██╔██╗██║██║░░██╗░█████╗░░██████╔╝  ██╔████╔██║██║░░██║██║░░██║█████╗░░██║░░░░░
░╚═══██╗██║░░██╗██╔══██║░░████╔═████║░██║██║╚████║██║░░╚██╗██╔══╝░░██╔══██╗  ██║╚██╔╝██║██║░░██║██║░░██║██╔══╝░░██║░░░░░
██████╔╝╚█████╔╝██║░░██║░░╚██╔╝░╚██╔╝░██║██║░╚███║╚██████╔╝███████╗██║░░██║  ██║░╚═╝░██║╚█████╔╝██████╔╝███████╗███████╗
╚═════╝░░╚════╝░╚═╝░░╚═╝░░░╚═╝░░░╚═╝░░╚═╝╚═╝░░╚══╝░╚═════╝░╚══════╝╚═╝░░╚═╝  ╚═╝░░░░░╚═╝░╚════╝░╚═════╝░╚══════╝╚══════╝
 * 
 * @author  Ambroise Muller
 * @date    November 2023
*/

#include "src/imports.hpp"


int main(int argc, char* argv[]) {

    if (argc == 1){

        // Random::rlxd_init(1, SEED);
        // Lattice* lattice = new Lattice(Nt, Nx, spacing);

        // PseudoFermionField *pf1 = lattice->pfermion[0];
        // PseudoFermionField *pf2 = lattice->pfermion[1];

        // lattice->gauge->values[0][0] = 0.3;
        // lattice->gauge->values[0][1] = 0.43;
        // lattice->gauge->values[1][0] = -10.3;
        // lattice->gauge->values[1][1] = 0.4;
        // lattice->gauge->values[2][0] = -7.4;
        // lattice->gauge->values[2][1] = 2.47;
        // lattice->gauge->values[3][0] = 3.1;
        // lattice->gauge->values[3][1] = -0.73;

        // pf1->values[0]->s[0] = complex<double>(0.2, -0.5);
        // pf1->values[0]->s[1] = complex<double>(-0.42, -0.25);
        // pf1->values[1]->s[0] = complex<double>(1.2, -0.45);
        // pf1->values[1]->s[1] = complex<double>(5.42, -0.95);
        // pf1->values[2]->s[0] = complex<double>(1.3, 0.6);
        // pf1->values[2]->s[1] = complex<double>(0.22, -0.23);
        // pf1->values[3]->s[0] = complex<double>(-0.14, -0.55);
        // pf1->values[3]->s[1] = complex<double>(0.16, 0.75);

        // ForceField *force = new ForceField(Nt, Nx);

        // pf1->compute_force2_to(force, 0.5, 0.1, 0.23, 1.0, 1e-10);
        
        // // pf2->set_to_dirac_of(pf1, 0.26, 0.5);
        // pf2->set_to_CG_solution(pf1, 0.26, 0.5, 1e-10, 1e4);
        // // pf2->add_z_times_pf(1.0, pf1);

        // for (int i=0; i<Nt*Nx; i++){
        //     cout << "site " << i << " : " << lattice->momentum->values[i][0] << ", " << lattice->momentum->values[i][1] << endl;
        // }

        // lattice->move_momentum(1.0);

        // for (int i=0; i<Nt*Nx; i++){
        //     cout << "site " << i << " : " << lattice->momentum->values[i][0] << ", " << lattice->momentum->values[i][1] << endl;
        // }

        // for (int i=0; i<Nt*Nx; i++){
        //     cout << "site " << i << endl;
        //     cout << pf1->values[i]->s[0].real() << " " << pf1->values[i]->s[0].imag() << " " << pf1->values[i]->s[1].real() << " " << pf1->values[i]->s[1].imag() << endl;
        // //     cout << pf2->values[i]->s[0].real() << " " << pf2->values[i]->s[0].imag() << " " << pf2->values[i]->s[1].real() << " " << pf2->values[i]->s[1].imag() << endl;
        //     cout << "force " << force->values[i][0] << ", " << force->values[i][1] << endl;
        // }


        // HoppingField *hop = lattice->gauge->hop;

        // cout << "gauge action: " << lattice->gauge->compute_gauge_action(4.0) << " qtop: " << lattice->gauge->compute_topological_charge() << endl;

        // lattice->gauge->compute_force_to(force, 4.0, 1.0);
        // for (int i=0; i<Nt*Nx; i++){
        //     cout << lattice->gauge->values[i][0] << " " << lattice->gauge->values[i][1] << " : " << force->values[i][0] << " " << force->values[i][1] << endl;
        // }
    
        Random::rlxd_init(1, SEED);
        Script::activate_environment();
        Lattice* lattice = new Lattice(Nt, Nx, spacing);
        lattice->HMC();
        lattice->gauge->save_config(OUTPUT_CONFIG);
        lattice->save_results(RESULTS_PATH);
        Script::run<Type::Python>("data_analysis/plot.py");
        Log::print(string("Results saved to folder ")+(RESULTS_PATH), Log::GOOD_NEWS);

    } else {
        Log::print("Number of arguments not supported!", Log::ERROR);
        exit(1);
    }

    return 0;
    
}
