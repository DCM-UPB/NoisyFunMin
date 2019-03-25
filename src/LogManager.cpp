#include "nfm/LogManager.hpp"

#include <fstream>
#include <iostream>
#include <sstream>


namespace nfm
{

LogLevel LogManager::log_level = LogLevel::OFF;
std::string LogManager::log_file_path;


void LogManager::setLoggingOn(const bool verbose)
{
    log_level = verbose ? LogLevel::VERBOSE : LogLevel::ON;
}

void LogManager::setLoggingOff()
{
    log_level = LogLevel::OFF;
}

void LogManager::setLogLevel(const LogLevel logLvl)
{
    log_level = logLvl;
}

bool LogManager::isLoggingOn()
{
    return (log_level > LogLevel::OFF);
}

bool LogManager::isVerbose()
{
    return (log_level > LogLevel::ON);
}

LogLevel LogManager::getLogLevel()
{
    return log_level;
}

void LogManager::setLoggingPathFile(const std::string &path)
{
    log_file_path = path;
}

void LogManager::writeOnLog(std::string s, const LogLevel logLvl)
{
    using namespace std;

    if (log_level < logLvl) { return; }

    const string log_marker = "--NFM--    ";
    const string log_marker_with_linebreak = "\n--NFM--    ";
    const string eol = "\n";

    if (log_file_path.empty()) {
        // append the string log_marker at the beginning of every new line
        size_t pos = s.find(eol, 0);
        // --- iterate through the string and change it accordingly
        while (pos < string::npos) {
            if (s.substr(pos, string::npos).length() > 0) {
                s.replace(pos, eol.length(), log_marker_with_linebreak);
            }
            pos = s.find(eol, pos + 1);
        }
        // write to the std output
        cout << log_marker << s << endl;
    }
    else {
        ofstream out;
        out.open(log_file_path, ios_base::app);
        out << s << endl;
        out.close();
    }
}

// --- Common logging routines

void LogManager::writeNoisyValueInLog(const NoisyValue nv, const LogLevel logLvl,
                                      const std::string &name, const std::string &flabel)
{
    using namespace std;

    if (log_level < logLvl) { return; }

    ostringstream os;
    if (!name.empty()) { os << name << ":\n"; }
    os << flabel << " = " << nv << "\n";
    os << flush;
    LogManager::writeOnLog(os.str(), logLvl);
}

void LogManager::writeVectorInLog(const std::vector<double> &x, const LogLevel logLvl,
                                  const std::string &name, const std::string &xlabel)
{
    using namespace std;

    if (log_level < logLvl) { return; }

    ostringstream os;
    if (!name.empty()) { os << name << ":\n"; }
    for (size_t i = 0; i < x.size(); ++i) {
        os << xlabel << i << " = " << x[i];
        if (i >= x.size() - 1) { os << "\n"; }
        else { os << "    "; }
    }
    os << flush;
    LogManager::writeOnLog(os.str(), logLvl);
}


void LogManager::writeNoisyVectorInLog(const std::vector<NoisyValue> &g, const LogLevel logLvl,
                                       const std::string &name, const std::string &glabel)
{
    using namespace std;

    if (log_level < logLvl) { return; }

    ostringstream os;
    if (!name.empty()) { os << name << ":\n"; }
    for (size_t i = 0; i < g.size(); ++i) {
        os << glabel << i << " = " << g[i];
        if (i >= g.size() - 1) { os << "\n"; }
        else { os << "    "; }
    }
    os << flush;
    LogManager::writeOnLog(os.str(), logLvl);
}


void LogManager::writeXValuePairInLog(const std::pair<std::vector<double>, NoisyValue> &pair, LogLevel logLvl,
                                      const std::string &name, const std::pair<std::string, std::string> &labels)
{
    LogManager::writeVectorInLog(pair.first, logLvl, name, labels.first);
    LogManager::writeNoisyValueInLog(pair.second, logLvl, "", labels.second);
}
} // namespace nfm