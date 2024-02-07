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
    void set_to(GaugeField *source){
        #if MULTI_TD
            auto set_to_lambda = [this, source](int start, int end, int thread_id) {
                for (int i = start; i < end; ++i) {
                    values[i][0] = source->values[i][0];
                    values[i][1] = source->values[i][1];
                }
            };
            HPC::MultithreadedLoop(set_to_lambda, 0, V);
        #else
            for (int i=0; i<V; i++){
                values[i][0] = source->values[i][0];
                values[i][1] = source->values[i][1];
            }
        #endif
    }

    /**
     * @brief Measure Wilson lattice action for current gauge field configuration.
     * @param beta Lattice coupling constant.
     * @return Contribution of gauge field to the total action.
    */
    double compute_gauge_action(double beta){
        #if MULTI_TD
            std::vector<double> partial_actions(std::thread::hardware_concurrency(), 0.0);
            auto action_lambda = [this, &partial_actions](int start, int end, int thread_id) {
                double local_action = 0.0;
                int imu, inu;
                double phi;
                for (int i = start; i < end; ++i) {
                    imu = hop->values[i][0];
                    inu = hop->values[i][1];
                    phi = values[i][0] + values[imu][1] - values[inu][0] - values[i][1];
                    local_action += cos(phi);
                }
                partial_actions[thread_id] = local_action;
            };
            HPC::MultithreadedLoop(action_lambda, 0, V);
            double total_action = std::accumulate(partial_actions.begin(), partial_actions.end(), 0.0);
            return -beta * total_action;
        #else
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
        #endif
    }

    /**
     * @brief Compute topological charge for current gauge field configuration.
     * @return Topological charge.
    */
    double compute_topological_charge(){
        #if MULTI_TD
            std::vector<double> partial_charges(std::thread::hardware_concurrency(), 0.0);
            double qtop_factor = 8.0 * atan(1.0);
            auto charge_lambda = [this, &partial_charges](int start, int end, int thread_id) {
                double local_qtop = 0.0;
                int imu, inu;
                double phi, offset;
                double qtop_factor = 8.0 * atan(1.0);
                for (int i = start; i < end; ++i) {
                    imu = (*hop)[i][0];
                    inu = (*hop)[i][1];
                    phi = values[i][0] + values[imu][1] - values[inu][0] - values[i][1];
                    offset = floor((phi + 0.5 * qtop_factor) / qtop_factor) * qtop_factor;
                    local_qtop += phi - offset;
                }
                partial_charges[thread_id] = local_qtop;
            };
            HPC::MultithreadedLoop(charge_lambda, 0, V);
            double total_qtop = std::accumulate(partial_charges.begin(), partial_charges.end(), 0.0);
            return total_qtop / qtop_factor;
        #else
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
        #endif
    }

    /**
     * @brief Evolve gauge field according to given momentum.
     * @param eps Step size.
     * @param mom Pointer to momentum field.
    */
    void move_gauge(double eps, MomentumField *mom) {
        #if MULTI_TD
            auto move_gauge_lambda = [this, eps, mom](int start, int end, int thread_id) {
                for (int i = start; i < end; ++i) {
                    for (int mu = 0; mu < D; ++mu) {
                        values[i][mu] -= eps * mom->values[i][mu];
                    }
                }
            };
            HPC::MultithreadedLoop(move_gauge_lambda, 0, V);
        #else 
            for (int i=0; i<V; i++){
                for (int mu=0; mu<D; mu++){
                    values[i][mu] -= eps * mom->values[i][mu];
                }
            }
        #endif
    }


    /**
     * @brief Compute gauge contribution to HMC force field.
     * @param force Pointer to ForceField instance where values should be written.
     * @param beta Lattice coupling constant.
     * @param a Constant @todo lattice spacing?
    */
    void compute_force_to(ForceField *force, double beta, double a){
        #if MULTI_TD
            auto compute_force_lambda = [this, force, beta, a](int start, int end, int thread_id) {
                int ix1, ix2;
                for (int i = start; i < end; ++i) {
                    for (int mu = 0; mu < D; mu++) {
                        for (int nu = 0; nu < D; nu++) {
                            if (mu == nu) continue;
                            // forward
                            ix1 = hop->values[i][mu];
                            ix2 = hop->values[i][nu];
                            force->values[i][mu] = -a * beta * sin(values[i][mu] + values[ix1][nu] - values[ix2][mu] - values[i][nu]);
                            // backward
                            ix1 = hop->values[i][D+nu];     
                            ix2 = hop->values[ix1][mu];  
                            force->values[i][mu] += -a * beta * sin(values[i][mu] - values[ix2][nu] - values[ix1][mu] + values[ix1][nu]);
                        }
                    }
                }
            };
            HPC::MultithreadedLoop(compute_force_lambda, 0, V);
        #else
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
        #endif
    }

};