
template <typename Type>
class Field {

public:

    int T;
    int L;
    int V;

    Type *values;

    Field(int T_, int L_) : T(T_), L(L_), V(T_*L_) {
        values = new Type[V];
    }

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


