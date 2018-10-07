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

    // --- Common logging routines

    void writeNoisyValueInLog(NoisyFunctionValue * x, const std::string &name = "Current position and target value", const std::string &xlabel = "x", const std::string &flabel = "f");

    void writeDirectionInLog(const double * grad, const int ndim, const double * dgrad = NULL /* optional errors */, const std::string &name = "Raw negative gradient", const std::string &glabel = "g");

    void reportMeaninglessGradientInLog();
};



#endif
