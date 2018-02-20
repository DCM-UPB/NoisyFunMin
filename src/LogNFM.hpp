#ifndef LOG_NFM
#define LOG_NFM

#include <string>


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

    void setLoggingPathFile(std::string path);

    void writeOnLog(std::string s);

};



#endif
