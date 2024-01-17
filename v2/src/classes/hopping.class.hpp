/**
 * @file    hopping.class.hpp
 * @brief   Nearest-neighbor field.
 *          Lattice point corresponding to coordinates (t, x) are at index [t*L + x].
 *          - HoppingField::values[t*L+x][0] = neighbour at (t, x+1)
 *          - HoppingField::values[t*L+x][1] = neighbour at (t+1, x)
 *          - HoppingField::values[t*L+x][2] = neighbour at (t, x-1)
 *          - HoppingField::values[t*L+x][3] = neighbour at (t-1, x)
 * 
 * @author  Ambroise Muller
 * @date    November 2023
*/


class HoppingField : public Field<int*> {

public:

    /**
     * @brief Constructor for HoppingField class.
     * @param T_ Temporal extent of lattice.
     * @param L_ Spatial extent of lattice. 
    */
    HoppingField(int T_, int L_) : Field(T_, L_) {
        for (int i=0; i<V; i++){
            values[i] = new int[2*D];
        }
    }

    /**
     * @brief Destructor for HoppingField class.
    */
    ~HoppingField(){
        for (int i=0; i<V; i++){
            delete values[i];
        }
    }

    /**
     * @brief Initialize neighbour field with periodic boundary conditions.
    */
    void initialize() {
        for (int i=0; i<V; i++){
            int x = i % L;
            int t = i / L;
            values[i][0] = t*L + ((x+1) % L);
            values[i][1] = ((t+1) % T)*L + x;
            values[i][2] = t*L + ((x-1+L) % L);
            values[i][3] = ((t-1+T) % T)*L + x;

        }
    }
};
