
#include "../structs/gauge.struct.hpp"

class GaugeField {

public:

    GaugeVector     **arr;
    int             T;
    int             L;
    int             V;

    GaugeField(int T_, int L_) : T(T_), L(L_), V(T_*L_) {
        arr = new GaugeVector*[V];
        for (int i=0; i<V; i++){
            arr[i] = new GaugeVector();
        }
    }

    ~GaugeField(){
        for (int i=0; i<V; i++){
            delete arr[i];
        }
        delete[] arr;
    }

    GaugeVector* operator[](int index){
        return arr[index];
    }

    const GaugeVector* operator[](int index) const {
        return arr[index];
    }

    void initialize_cold(){
        for (int i=0; i<V; i++){
            arr[i]->t = 0;
            arr[i]->x = 0;
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
            arr[i]->t = std::stod(token);
            std::getline(iss, token, ',');
            arr[i]->x = std::stod(token);
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
            output << to_string(arr[i]->t) << ", " << to_string(arr[i]->x) << endl;
        }
        Log::print("Successfully saved gauge field configuration to: "+filename, Log::GOOD_NEWS);
        output.close();
    }

};