#include "nfm/LogNFM.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <thread>

namespace nfm
{

int NFMLogManager::_log_level = 0;
std::string NFMLogManager::_log_file_path;


void NFMLogManager::setLoggingOn(const bool verbose)
{
    _log_level = verbose ? 2 : 1;
}

void NFMLogManager::setLoggingOff()
{
    _log_level = 0;
}

void NFMLogManager::setLogLevel(const int log_level)
{
    _log_level = log_level;
}

bool NFMLogManager::isLoggingOn()
{
    return (_log_level > 0);
}

bool NFMLogManager::isVerbose()
{
    return (_log_level > 1);
}

int NFMLogManager::getLogLevel()
{
    return _log_level;
}

void NFMLogManager::setLoggingPathFile(const std::string &path)
{
    _log_file_path=path;
}

void NFMLogManager::writeOnLog(std::string s, const int log_level)
{
    using namespace std;

    if (_log_level < log_level) { return; }

    const std::string log_marker = "--NFM--    ";
    const std::string log_marker_with_linebreak = "\n--NFM--    ";

    if (_log_file_path.empty()){
        // append the string log_marker at the beginning of every new line
        size_t pos = s.find('\n', 0);
        // --- iterate through the string and change it accordingly
        while(pos < string::npos){
            if (s.substr(pos, string::npos).length() > 0) {
                s.replace(pos, string("\n").length(), log_marker_with_linebreak);
            }
            pos = s.find('\n', pos+1);
        }
        // write to the std output
        cout << log_marker << s << endl;
    } else {
        ofstream out;
        out.open(_log_file_path, ios_base::app);
        out << s << endl;
        out.close();
    }

    //this_thread::sleep_for (chrono::milliseconds(150));
}

// --- Common logging routines

void NFMLogManager::writeNoisyValueInLog(NoisyFunctionValue * x, const int log_level, const std::string &name, const std::string &flabel, const bool printX, const std::string &xlabel)
{
    using namespace std;

    if (_log_level<log_level) { return; }

    stringstream s;
    s << name << ":\n";
    if (printX) {
        for (int i=0; i<x->getNDim(); ++i){
            s << xlabel << i << " = " << x->getX(i);
            if (i==x->getNDim()-1) { s << endl; }
            else { s << "    "; }
        }
    }
    s << flabel << " = " << x->getF() << " +- " << x->getDf() << endl;
    s << flush;
    this->writeOnLog(s.str(), log_level);
}


void NFMLogManager::writeVectorInLog(const double * vec, const double * dvec, const int ndim, const int log_level, const std::string &name, const std::string &vlabel)
{
    using namespace std;

    if (_log_level<log_level) { return; }

    stringstream s;
    s << name << ":\n";
    for (int i=0; i<ndim; ++i){
        s << vlabel << i << " = " << vec[i];
        if (dvec!=nullptr) { s << " +- " << dvec[i]; }
        if (i==ndim-1 ) { s << endl; }
        else { s << "    "; }
    }
    s << flush;
    this->writeOnLog(s.str(), log_level);
}
} // namespace nfm