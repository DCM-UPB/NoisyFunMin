#ifndef NFM_NOISYFUNMIN_HPP
#define NFM_NOISYFUNMIN_HPP

#include <list>
#include <string>

#include "nfm/NoisyFunction.hpp"
#include "nfm/NoisyValue.hpp"

namespace nfm
{

class NFM
{
protected:
    const int _ndim;  //dimensionality of the space where the target function is embedded
    NoisyFunction * const _targetfun;  //target function to minimize

    NoisyFunctionWithGradient * const _gradfun;  //gradient of the target function
    const bool _flag_gradfun;  //has the gradient been provided?
    const bool _flag_graderr; // use a provided gradient error for printout / stopping criteria

    const int _max_n_const_values; //stop after this number of target values have been constant within error bounds
    std::list<NoisyValue> _old_values; // list of previous target values

    NoisyIOPair _last;  //last position and its corresponding function value

    double _epsf; //changes in the function smaller than this value will stop the minimization
    double _epsx; //changes in the position x smaller than this value will stop the minimization

    void _clearOldValues() { _old_values.clear(); } // reset old values list
    void _storeOldValue(); // store last value in old values list
    bool _isConverged() const; // check if the target function has stabilized
    bool _meaningfulGradient(const std::vector<NoisyValue> &grad) const; //check if the gradient is meaningful. i.e. if its values are greater than the statistical errors
    bool _shouldStop(const std::vector<NoisyValue> * grad) const; // check for all stopping criteria

    // "Mandatory" logging routines
    // CurrentX and XUpdate should be logged by all optimizers and also GradientInLog if a gradient is used
    void _writeCurrentXToLog() const;
    void _writeGradientToLog(const std::vector<NoisyValue> &grad) const;
    void _writeXUpdateToLog(const std::vector<double> &xu) const;

    // Stopping criterium debug logger
    void _writeOldValuesToLog();

public:
    explicit NFM(NoisyFunction * targetfun, int max_n_const_values = 20);
    virtual ~NFM() = default;

    // --- Setters
    void setX(int i, double x) { _last.x[i] = x; } // set per element
    void setX(const double x[]); // set via c-style array
    void setX(const std::vector<double> &x) { this->setX(x.data()); } // set via vector
    void setEpsF(double epsf) { _epsf = epsf; }
    void setEpsX(double epsx) { _epsx = epsx; }

    // --- Getters

    int getNDim() const { return _ndim; }
    NoisyFunction * getTargetFun() const { return _targetfun; }
    NoisyFunctionWithGradient * getGradientFun() const { return _gradfun; }

    double getX(int i) const { return _last.x[i]; }; // elementary get
    void getX(double x[]) const; // get via c-style array
    void getX(std::vector<double> &x) const { return this->getX(x.data()); } // get via vector
    std::vector<double> getX() const { return _last.x; } // get in new vector
    double getF() const { return _last.f.value; }
    double getDf() const { return _last.f.error; }
    NoisyValue getFDf() const { return _last.f; }

    double getEpsF() const { return _epsf; }
    double getEpsX() const { return _epsx; }

    // --- Minimization
    virtual void findMin() = 0;
};
} // namespace nfm

#endif
