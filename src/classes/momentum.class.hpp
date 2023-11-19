/**
 * @file    momentum.class.hpp
 * @brief   Momentum field for HMC simulations of the lattice Schwinger model.
 *          Lattice point corresponding to coordinates (t, x) are at index [t*L + x].
 * 
 * @author  Ambroise Muller
 * @date    November 2023
*/

class MomentumField : public Field<double*> {

public:

    MomentumField(int T_, int L_) : Field(T_, L_) {
        for (int i=0; i<V; i++){
            values[i] = new double[D];
        }
    }

    ~MomentumField() {
        for (int i=0; i<V; i++){
            delete values[i];
        }
    }

    double initialize_momentum(){
        double norm = 0.0;
        for (int i=0; i<V; i++){
            for (int mu=0; mu<D; mu++){
                Random::gauss_rand(1, values[1]+mu);
                values[i][mu] *= sqrt(2);
                norm += values[i][mu] * values[i][mu];
            }
        }
        return 0.5*norm;
    }

    double compute_momentum_action(){
        double action = 0.0;
        for (int i=0; i<V; i++){
            for (int mu=0; mu<D; mu++){
                action += values[i][mu] * values[i][mu];
            }
        }
        return 0.5*action;
    }

};