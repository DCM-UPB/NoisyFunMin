#ifndef LOG_NFM
#define LOG_NFM

#include <string>

#include "nfm/NoisyFunctionValue.hpp"


class NFMLogManager{
private:
    static int _log_level; // <1 = no log, 1 = essential log, >1 verbose log
    static std::string _log_file_path;  // if the path is not given, the log uses the cout

public:
    NFMLogManager(){}
    ~NFMLogManager(){}

    void setLoggingOn(const bool verbose = false); // set loglevel to 1 or 2 (verbose)
    void setLoggingOff(); // set loglevel 0
    void setLogLevel(const int log_level);
    bool isLoggingOn();
    bool isVerbose(); // log_level > 1?
    int getLogLevel();

    void setLoggingPathFile(const std::string &path);

    void writeOnLog(std::string s, const int log_level = 1);

    // --- Advanced log helpers
    void writeNoisyValueInLog(NoisyFunctionValue * x, const int log_level, const std::string &name, const std::string &flabel = "f", const bool printX = false, const std::string &xlabel = "x");
    void writeVectorInLog(const double * vec, const double * dvec /* optional error vector */, const int ndim, const int log_level, const std::string &name, const std::string &vlabel = "v");
};



#endif
