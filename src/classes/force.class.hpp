/**
 * @file    force.class.hpp
 * @brief   Force fields used for molecular dynamics in HMC.
 *          Lattice point corresponding to coordinates (t, x) are at index [t*L + x].
 * 
 * @author  Ambroise Muller
 * @date    November 2023
*/

class ForceField : public Field<double*> {

public:

    /**
     * @brief Constructor for ForceField class.
     * @param T_ Temporal extent of lattice.
     * @param L_ Spatial extent of lattice. 
    */
    ForceField(int T_, int L_) : Field(T_, L_) {
        for (int i=0; i<V; i++){
            values[i] = new double[D];
        }
    }

    /**
     * @brief Destructor for ForceField class.
    */
    ~ForceField() {
        for (int i=0; i<V; i++){
            delete values[i];
        }
    }

    /**
     * @brief Set all forces to zero.
    */
    void set_zero(){
        for (int i=0; i<V; i++){
            for (int mu=0; mu<D; mu++){
                values[i][mu] = 0.0;
            }
        }
    }

    /**
     * @brief Overload of in-place addition.
     * @param other Pointer to other ForceField instance.
     * @return Pointer to current instance.
    */
    ForceField& operator+=(const ForceField& other) {
        if (this == &other) {
            return *this;
        }
        for (int i = 0; i < V; i++) {
            for (int mu = 0; mu < D; mu++) {
                this->values[i][mu] += other.values[i][mu];
            }
        }
        return *this; 
    }

};