#include "LogNFM.hpp"

#include <iostream>
#include <fstream>
#include <thread>



bool NFMLogManager::_is_logging_on = false;
std::string NFMLogManager::_log_file_path = "";



void NFMLogManager::setLoggingOn(){
    _is_logging_on=true;
}


void NFMLogManager::setLoggingOff(){
    _is_logging_on=false;
}


bool NFMLogManager::isLoggingOn(){
    return _is_logging_on;
}


void NFMLogManager::setLoggingPathFile(std::string path){
    _log_file_path=path;
}


void NFMLogManager::writeOnLog(std::string s){
    using namespace std;

    if (this->isLoggingOn()){
        const std::string log_marker = "--NFM--    ";
        const std::string log_marker_with_linebreak = "\n--NFM--    ";

        if (_log_file_path.empty()){
            // append the string log_marker at the beginning of every new line
            size_t pos = s.find("\n", 0);
            // --- iterate through the string and change it accordingly
            while(pos < string::npos){
                if (s.substr(pos, string::npos).length() > 0)
                    s.replace(pos, string("\n").length(), log_marker_with_linebreak);
                pos = s.find("\n", pos+1);
            }
            // write to the std output
            cout << log_marker << s << endl;
        } else{
            ofstream out;
            out.open(_log_file_path, ios_base::app);
            out << s << endl;
            out.close();
        }

        //this_thread::sleep_for (chrono::milliseconds(150));
    }
}
