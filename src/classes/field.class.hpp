/**
 * @file    field.class.hpp
 * @brief   Template class for fields defined on the lattice.
 *          Lattice point corresponding to coordinates (t, x) are at index [t*L + x].
 * 
 * @author  Ambroise Muller
 * @date    November 2023
*/


template <typename Type>
class Field {

public:

    const int   T;          // Number of sites in temporal direction
    const int   L;          // Number of sites in spatial direction
    const int   V;          // Lattice volume

    Type *values;           // Field values

    /**
     * @brief Constructor for Field class.
     * @param T_ Temporal extent of lattice.
     * @param L_ Spatial extent of lattice. 
    */
    Field(int T_, int L_) : T(T_), L(L_), V(T_*L_) {
        values = new Type[V];
    }

    /**
     * @brief Destructor for Field class.
    */
    ~Field(){
        delete[] values;
    }

    /**
     * @brief Overload of subscript operator.
     * @param index Index of lattice point to fetch.
    */
    Type& operator[](int index){
        return values[index];
    }

    /**
     * @brief Read-only overload of subscript operator.
     * @param index Index of lattice point to fetch.
    */
    const Type& operator[](int index) const {
        return values[index];
    }

};


