/**
 * @file    gauge.class.hpp
 * @brief   Gauge field for the Schwinger model.
 *          Lattice point corresponding to coordinates (t, x) are at index [t*L + x].
 * 
 * @author  Ambroise Muller
 * @date    November 2023
*/

class GaugeField : public Field<double*> {

public:

    HoppingField    *hop;   // Hopping field

    /**
     * @brief Constructor for GaugeField class.
     * @param T_ Temporal extent of lattice.
     * @param L_ Spatial extent of lattice. 
    */
    GaugeField(int T_, int L_) : Field(T_, L_) {
        for (int i=0; i<V; i++){
            values[i] = new double[D];
        }
    }

    /**
     * @brief Destructor for GaugeField class.
    */
    ~GaugeField(){
        for (int i=0; i<V; i++){
            delete values[i];
        }
    }

    /**
     * @brief Assign hopping field.
     * @param hopping Pointer to hopping field.
    */
    void assign_hopping_field(HoppingField *hopping){
        hop = hopping;
        Log::print("Successfully assigned hopping field", Log::VERBOSE);
    }

    /**
     * @brief Initialize the gauge field in a cold (i.e. trivial) configuration.
    */
    void initialize_cold(){
        for (int i=0; i<V; i++){
            values[i][0] = 0;
            values[i][1] = 0;
        }
    }

    /**
     * @brief Load gauge field configuration for file.
     * @param filename Path to file holding saved gauge field configuration.
    */
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

    /**
     * @brief Save current gauge field configuration to file.
     * @param filename Path to file to hold saved gauge field configuration.
    */
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

    /**
     * @brief Assign the values of a source field to current gauge field.
     * @param source Pointer to source gauge field.
    */
    void assign_link(GaugeField *source){
        for (int i=0; i<V; i++){
            values[i][0] = source->values[i][0];
            values[i][1] = source->values[i][1];
        }
    }

    /**
     * @brief Measure Wilson lattice action for current gauge field configuration.
     * @param beta Lattice coupling constant.
     * @return Contribution of gauge field to the total action.
    */
    double compute_gauge_action(double beta){
        int imu, inu;
        double phi;
        double action = 0.0;
        for (int i=0; i<V; i++){
            imu = hop->values[i][0];
            inu = hop->values[i][1];
            phi = values[i][0] + values[imu][1] - values[inu][0] - values[i][1];
            action += cos(phi);
        }
        return -beta * action;
    }

    /**
     * @brief Compute topological charge for current gauge field configuration.
     * @return Topological charge.
    */
    double compute_topological_charge(){
        int imu, inu;
        double phi, offset;
        double qtop = 0.0;
        double qtop_factor = 8.0*atan(1.0);
        for (int i=0; i<V; i++){
            imu = (*hop)[i][0];
            inu = (*hop)[i][1];
            phi = values[i][0] + values[imu][1] - values[inu][0] - values[i][1];
            offset = floor((phi+0.5*qtop_factor)/qtop_factor) * qtop_factor;
            qtop += phi - offset;
        }
        return qtop/qtop_factor;
    }

    /**
     * @brief Evolve gauge field according to given momentum.
     * @param eps Step size.
     * @param mom Pointer to momentum field.
    */
    void move_gauge(double eps, MomentumField *mom){
        for (int i=0; i<V; i++){
            for (int mu=0; mu<D; mu++){
                values[i][mu] -= eps * mom->values[i][mu];
            }
        }
    }

    /**
     * @brief Compute gauge contribution to HMC force field.
     * @param force Pointer to ForceField instance where values should be written.
     * @param beta Lattice coupling constant.
     * @param a Constant @todo lattice spacing?
    */
    void compute_force_to(ForceField *force, double beta, double a){
        int ix1, ix2;
        for (int i=0; i<V; i++){
            for (int mu=0; mu<D; mu++){
                for (int nu=0; nu<D; nu++){
                    if (mu == nu) {continue;}
                    // forward
                    ix1 = hop->values[i][mu];
                    ix2 = hop->values[i][nu];
                    force->values[i][mu] = -a * beta * sin(values[i][mu] + values[ix1][nu] - values[ix2][mu] - values[i][nu]);
                    // backward
                    ix1 = hop->values[i][D+nu];     
                    ix2 = hop->values[ix1][mu];     // @todo this looks fishy (ix1->i, mu->D+mu ?)
                    force->values[i][mu] += -a * beta * sin(values[i][mu] - values[ix2][nu] - values[ix1][mu] + values[ix1][nu]);
                }
            }
        }
    }

};