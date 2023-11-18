

class GaugeField {

public:

    double  **values;   // stores gauge field values for each lattice site
    int     T;          // number of sites along temporal direction
    int     L;          // number of sites along spatial direction
    int     V;          // total lattice volume

    GaugeField(int T_, int L_) : T(T_), L(L_), V(T_*L_) {
        values = new double*[V];
        for (int i=0; i<V; i++){
            values[i] = new double[D];
        }
    }

    ~GaugeField(){
        for (int i=0; i<V; i++){
            delete values[i];
        }
        delete[] values;
    }

    double* operator[](int index){
        return values[index];
    }

    const double* operator[](int index) const {
        return values[index];
    }

    void initialize_cold(){
        for (int i=0; i<V; i++){
            values[i][0] = 0;
            values[i][0] = 0;
        }
    }

    void initialize_from(const string& filename){
        ifstream input(filename);
        if (!input.is_open()) {
            Log::print("Error: unable to open file for reading at path: "+filename, Log::ERROR);
            exit(1);
        }
        string line;
        for (int i = 0; i < V && std::getline(input, line); i++) {
            istringstream iss(line);
            string token;
            std::getline(iss, token, ',');
            values[i][0] = std::stod(token);
            std::getline(iss, token, ',');
            values[i][1] = std::stod(token);
        }
        input.close();
    }

    void save_config(const string& filename){
        ofstream output(filename);
        if (!output.is_open()) {
            Log::print("Error: Unable to open file for writing at path: "+filename, Log::ERROR);
            exit(1);
        }
        for (int i=0; i<V; i++){
            output << to_string(values[i][0]) << ", " << to_string(values[i][1]) << endl;
        }
        Log::print("Successfully saved gauge field configuration to: "+filename, Log::GOOD_NEWS);
        output.close();
    }

};