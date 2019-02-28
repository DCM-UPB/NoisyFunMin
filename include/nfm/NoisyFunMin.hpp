#ifndef NFM_NOISYFUNMIN_HPP
#define NFM_NOISYFUNMIN_HPP

#include <list>
#include <string>

#include "nfm/NoisyFunction.hpp"
#include "nfm/NoisyFunctionValue.hpp"


class NFM
{
protected:
    int _ndim;  //dimensionality of the space where the target function is embedded
    NoisyFunction * _targetfun;  //target function to minimize

    NoisyFunctionValue * _x;  //last position and its corresponding function value

    NoisyFunctionWithGradient * _gradtargetfun;  //gradient of the target function
    bool _flaggradtargetfun;  //has the gradient been provided?
    const bool _useGradientError; // use a provided gradient error for printout / stopping criteria

    //nfm::DomainFun _indomain; //function that determines if a point is inside the domain of optimization or not
    //bool _flagindomain; //has the domain been provided?

    double _epstargetfun; //changes in the function smaller than this value will stop the minimization
    double _epsx; //changes in the position x smaller than this value will stop the minimization

    const size_t _max_n_const_values; //stop after this number of target values have been constant within error bounds
    std::list<NoisyFunctionValue *> _old_values; // list of previous target values

    void _clearOldValues(); // reset old values list
    void _storeOldValue(); // store last value in old values list
    bool _isConverged(); // check if the target function has stabilized
    bool _meaningfulGradient(const double * grad, const double * graderr); //check if the gradient is meaningful. i.e. if its values are greater than the statistical errors
    bool _shouldStop(const double * grad, const double * graderr); // check for all stopping criteria

    // "Mandatory" logging routines
    // CurrentX and XUpdate should be logged by all optimizers and also GradientInLog if a gradient is used
    void _writeCurrentXInLog();
    void _writeGradientInLog(const double * grad, const double * dgrad);
    void _writeXUpdateInLog(const double * xu);

    // Stopping criterium debug logger
    void _writeOldValuesInLog();

public:
    explicit NFM(NoisyFunction * targetfun, bool useGradientError = true, const size_t &max_n_const_values = 20);
    virtual ~NFM();

    // --- Setters
    void setX(const double * x);
    void setX(const int &i, const double &x);
    void setGradientTargetFun(NoisyFunctionWithGradient * grad);
    //void setDomain(nfm::DomainFun domain);
    void setEpsTargetFun(double &epstargetfun){_epstargetfun=epstargetfun;}
    void setEpsX(double &epsx){_epsx=epsx;}

    // --- Getters
    int getNDim(){return _ndim;}
    double getX(const int &i){return _x->getX(i);};
    void getX(double * x){for(int i=0; i<_ndim; ++i){x[i]=_x->getX(i);}}
    double getF(){ return _x->getF(); }
    double getDf(){ return _x->getDf(); }
    NoisyFunctionWithGradient* getGradientTargetFun(){return _gradtargetfun;}
    //nfm::DomainFun getDomain(){return _indomain;}
    double getEpsTargetFun(){return _epstargetfun;}
    double getEpsX(){return _epsx;}

    // --- Minimization
    virtual void findMin() = 0;
};


#endif
