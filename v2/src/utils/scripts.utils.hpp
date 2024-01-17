/**
 * @file
 * @brief
 * 
 * @author  Ambroise Muller
 * @date    November 2023
*/

enum class Type {
    Python
};

/**
 * @brief The script namespace
*/
namespace Script {

    string data_folder = "";

    void activate_environment(){
        #ifdef _WIN32
            system("../venv/Scripts/activate");
        #else
            system("source ../venv/bin/activate");
        #endif
    }

    template<Type T>
    void run(string script_name, string blend_file);

    template<Type T>
    void run(string script_name);

    /**
     * @brief Run a python script
     * @param script_name The name of the script
    */
    template<>
    void run<Type::Python>(string script_name) {
        Log::title("Run python script");
        #ifdef _WIN32
            system(("../venv/Scripts/python ../" + script_name).c_str());
        #else
	        system(("../venv/bin/python ../" + script_name).c_str());
	    #endif
    }
};