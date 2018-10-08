#ifndef LOG_NFM
#define LOG_NFM

#include <string>

#include "NoisyFunctionValue.hpp"


class NFMLogManager{
private:
    static bool _is_logging_on;
    static std::string _log_file_path;  // if the path is not given, the log uses the cout

public:
    NFMLogManager(){}
    ~NFMLogManager(){}

    void setLoggingOn();
    void setLoggingOff();
    bool isLoggingOn();

    void setLoggingPathFile(const std::string &path);

    void writeOnLog(std::string s);

    // --- Advanced log helpers
    void writeNoisyValueInLog(NoisyFunctionValue * x, const std::string &name, const std::string &xlabel = "x", const std::string &flabel = "f");
    void writeVectorInLog(const double * vec, const double * dvec /* optional error vector */, const int ndim, const std::string &name, const std::string &vlabel = "v");
};



#endif
