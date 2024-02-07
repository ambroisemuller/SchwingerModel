/**
 * @file    log.utils.hpp
 * @brief   Define Log namespace and associated functionality.
 * 
 * @author  Ambroise Muller
 * @date    November 2023
*/

#include "../configs/log.config.h"

#define lineWidth 100

class ProgressBar {
    public:
        explicit ProgressBar() : total_(0), progress_(0), mutex_() {}

        void increment() {
            std::lock_guard<std::mutex> lock(mutex_);
            ++progress_;
            update();
        }

        void set_total(double total) {
            total_ = total;
        }

    private:
        int total_;
        int progress_;
        std::mutex mutex_;

        void update() {
            double percentage = static_cast<double>(progress_) / total_ * 100;
            int filledCount = static_cast<int>(percentage* 72./100.);
            int emptyCount = lineWidth - filledCount - 7;

            std::cout << "[";
            for (int i = 0; i < filledCount; ++i)
                std::cout << "=";
            for (int i = 0; i < emptyCount; ++i)
                std::cout << " ";
            std::cout << "] " << static_cast<int>(percentage) << "%\r";
            std::cout.flush();

            if (progress_ == total_)
                std::cout << std::endl;
        }
};

/**
 * @brief The Log namespace
*/
namespace Log {

    typedef int LOG_LEVEL_td;
    typedef std::string LOG_NAME_td;

    std::stringstream* log = new std::stringstream;

    const LOG_LEVEL_td GOOD_NEWS  = 1;   // Used for reporting good news, the highest priority of content
    const LOG_LEVEL_td ERROR      = 1;   // Used for reporting errors, the highest priority of content
    const LOG_LEVEL_td WARNING    = 2;   // Used for warnings. Device should still be partly working (and making do).
    const LOG_LEVEL_td VERBOSE    = 3;   // Used for dumping everything to print. Also check

    double t;                            // Time
    std::string bar_label = "";
    auto current_time = std::chrono::high_resolution_clock::now();

    bool activ = true;

    ProgressBar progress_bar = ProgressBar();

    /**
     * @brief Set the chrono to the current time to then find the difference
    */
    void init_time() {
        current_time = std::chrono::high_resolution_clock::now();
    }

    /**
     * @brief Print a progress bar
     * @param progress The state of the progress bar ranging from 0 to 1
     * @param label bool Print bar_label in front of progress bar
    */
    void progressBar(double progress, bool label=false) {

        int elapsed_time = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::high_resolution_clock::now() - current_time).count();

        if (label) {
            std::cout << "\r" << bar_label;

            for (int i = 0; i < 8 - (int)(bar_label.size()); ++i) {
                std::cout << " ";
            }
        }
        else {
            std::cout << "\r";
        }
        std::cout << "[";

        int pos = int(progress * (lineWidth-73 - (label? 8:0)));

        for (int i = 0; i < (lineWidth-73 - (label? 8:0)); ++i) {
            if (i <= pos) std::cout << "=";
            else std::cout << " ";
        }

        std::cout << "] " << ceil(progress * 100.0) << " % | Elapsed: " << std::to_string(int(elapsed_time)) << "s"
        << " Est. remaining: " << std::to_string(max(0, int(elapsed_time * (1 / (progress+1e-64) - 1)))) << "s"
        << " Est. total: " << std::to_string(max(0, int(elapsed_time / (progress+1e-64) ))) << "s";

        std::cout.flush();
    }

    /**
     * @brief Print a nice green title
    */
    void title(const std::string& str) {

        std::cout << std::endl << "\033[1;32m";

        for (int i = 0; i < int((lineWidth - str.length())/2 - 1); ++i) {
            std::cout << "-";
        }
        std::cout << " " << str << " ";

        for (int i = 0; i < int((lineWidth - str.length() + 1)/2 - 1); ++i) {
            std::cout << "-";
        }

        std::cout << "\033[0m" << std::endl;
    }

    #if DEBUG

        /**
         * @brief Print the value of a variable with a label
         * @param str : the label of the variable
         * @param variable : the value of the variable
        */
        template<typename T>
        void print(const std::string& str, const T& variable, LOG_LEVEL_td level = VERBOSE) {
            if (level <= LOG_LEVEL && activ) {
                std::cout << "[" << to_string(t) << "] " << str << ": " << get_string(variable) << std::endl;
            }
            if (SAVE_LOG && level <= SAVE_LEVEL && activ) {
                *log << "[" << to_string(t) << "] " << str << ": " << get_string(variable) << std::endl;
            }
        }

        /**
         * @brief Simple print
         * @param str : the text
        */
        void print(const std::string& str, LOG_LEVEL_td level = VERBOSE) {
            if (level <= LOG_LEVEL && activ) {
                std::cout << "[" << to_string(t) << "] " << str << std::endl;
            }
            if (SAVE_LOG && level <= SAVE_LEVEL && activ) {
                *log << "[" << to_string(t) << "] " << str << std::endl;
            }
        }

    #else

        template<typename T>
        void print(const std::string& str, const T& variable, LOG_LEVEL_td level = VERBOSE) {}
        void print(const std::string& str, LOG_LEVEL_td level = VERBOSE) {}

    #endif

    /**
     * @brief Go to the next line
    */
    void newLine() {
        if (activ) {
            std::cout << std::endl;
        }
        if (SAVE_LOG && activ) {
            *log << std::endl;
        }
    }

    /**
     * @brief Save the log with the specified name, or in '../results/'log.txt' by default
     * @param folder Optional name of log file
     */
    void save(const string& folder = "results/log.txt") {
        std::string filename = "../" + folder;
        std::ofstream output(filename);
        if (!output.is_open()) {
            std::cerr << "Error: Unable to open file for writing at path: " << filename << std::endl;
            return;
        }
        output << log->str();
    }
}