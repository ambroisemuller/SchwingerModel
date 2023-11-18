

class HoppingField {

public:

    int **values;   // stores gauge field values for each lattice site
    int T;          // number of sites along temporal direction
    int L;          // number of sites along spatial direction
    int V;          // total lattice volume

    HoppingField(int T_, int L_) : T(T_), L(L_), V(T_*L_) {
        values = new int*[V];
        for (int i=0; i<V; i++){
            values[i] = new int[2*D];
        }
    }

    ~HoppingField(){
        for (int i=0; i<V; i++){
            delete values[i];
        }
        delete[] values;
    }

    int* operator[](int index){
        return values[index];
    }

    const int* operator[](int index) const {
        return values[index];
    }

    /**
     * @brief Initialize neighbour field with periodic boundary conditions
    */
    void initialize() {
        cout << "here" << endl;
        cout << "here2";
        for (int x=0; x<V; x++){
            int Lk = V;
            int y = x;
            for (int k=D-1; k>=0; k--){
                Lk /= L;                       
                int xk = y/Lk;                    
                y  = y - xk * Lk;          
                int dxk;
                // forward
                if (xk < L-1){
                    dxk = Lk;
                } else {        
                    dxk = Lk * (1-L);
                }
                values[x][k] = x + dxk;
                // backward
                if (xk > 0){   
                    dxk = -Lk;
                } else {        
                    dxk = Lk * (L-1);
                }
                values[x][k+D] = x + dxk;
            }
        }
        cout << "there" << endl;
    }

};