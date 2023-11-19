

class ForceField : public Field<double*> {

public:

    ForceField(int T_, int L_) : Field(T_, L_) {
        for (int i=0; i<V; i++){
            values[i] = new double[D];
        }
    }

    ~ForceField() {
        for (int i=0; i<V; i++){
            delete values[i];
        }
    }

    void set_zero(){
        for (int i=0; i<V; i++){
            for (int mu=0; mu<D; mu++){
                values[i][mu] = 0.0;
            }
        }
    }

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