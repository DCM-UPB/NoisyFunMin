#ifndef NFM_LOGMANAGER_HPP
#define NFM_LOGMANAGER_HPP

#include "nfm/NoisyValue.hpp"

#include <string>
#include <vector>
#include <tuple>

namespace nfm
{

// settable log levels
enum class LogLevel
{
    OFF = 0,
    ON = 1,
    VERBOSE = 2
};

// For simplicity we currently use a LogManager with
// static fields. This means very easy output control
// in application code, without requiring direct access
// to a NFM object. However, it may be problematic in
// certain multi-threading scenarios with shared memory.
class LogManager
{
public:
    static LogLevel log_level; // <1 = no log, 1 = essential log, >1 verbose log
    static std::string log_file_path;  // if the path is not given, the log uses the cout

    static void setLoggingOn(bool verbose = false); // set loglevel to ON or VERBOSE
    static void setLogLevel(LogLevel level = LogLevel::ON); // passing ints should work too
    static void setLoggingOff();
    static LogLevel getLogLevel();
    static bool isLoggingOn();
    static bool isVerbose();

    static void setLoggingPathFile(const std::string &path);

    static void writeOnLog(std::string s, LogLevel logLvl = LogLevel::ON);

    // --- Advanced log helpers

    // write single noisy value (e.g. NoisyFunction result)
    static void writeNoisyValueInLog(NoisyValue nv, LogLevel logLvl,
                                     const std::string &name = "", const std::string &flabel = "f");

    // write vector of exact values (e.g. positions)
    static void writeVectorInLog(const std::vector<double> &x, LogLevel logLvl,
                                 const std::string &name = "", const std::string &xlabel = "x");

    // write vector of noisy values (e.g. gradients)
    static void writeNoisyVectorInLog(const std::vector<NoisyValue> &g, LogLevel logLvl,
                                      const std::string &name = "", const std::string &glabel = "g");

    // write a pair of exact vector and noisy value
    static void writeXValuePairInLog(const std::pair<std::vector<double>, NoisyValue> &pair, LogLevel logLvl,
                                     const std::string &name = "", const std::pair<std::string, std::string> &labels = {"x", "f"});
};
} // namespace nfm

#endif
