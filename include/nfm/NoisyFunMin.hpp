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
    // Default values for general NFM parameters
    static constexpr int DEFAULT_MAX_N_CONST = 20; // stop on maximal number of constant values
    static constexpr double DEFAULT_STEPSIZE = 0.01; // default to use for step sizes
    static constexpr double DEFAULT_EPSX = 1.e-5; // stop on minimal position change
    static constexpr double DEFAULT_EPSF = 0.; // stop on minimal target function change

    // Protected members
    const int _ndim;  //dimensionality of the space where the target function is embedded
    NoisyFunction * const _targetfun;  //target function to minimize

    NoisyFunctionWithGradient * const _gradfun;  //gradient of the target function
    const bool _flag_gradfun; // has the gradient been provided?
    const bool _flag_graderr; // does the gradfun provide gradient errors?

    // Tolerance (may also be used by child to auto-set own tolerances)
    double _epsx; // changes in the position x smaller than this value will stop the minimization
    double _epsf; // changes in the function smaller than this value will stop the minimization

    NoisyIOPair _last;  // last position and its function value, to be updated by child

private: // private methods used by base class
    const int _max_n_const_values; // stop after this number of target values have been constant within error bounds
    std::list<NoisyIOPair> _old_values; // list of previous target values
    double _lastDeltaX{}; // change in x by last step (updated on storeLastValue)
    double _lastDeltaF{}; // change in f by last step (updated on storeLastalue)

    void _clearOldValues() { _old_values.clear(); } // reset old values list
    bool _isConverged() const; // check if the target function has stabilized
    void _updateDeltas(); // calculate deltaX and deltaF between _last and _old_values.front()
    bool _checkDeltas() const; // check deltas against epsx and epsf

    void _writeCurrentXToLog() const; // write current x on log on storeLastValue
    void _writeOldValuesToLog() const; // stopping criterium debug logger

protected: // protected methods for derived optimizers
    // use this after every position&function update
    void _storeLastValue(); // store last value in old values list (updates deltax/deltaf)

    // check stopping criteria (shouldStop contains meaningfulGradient, if grad!=nullptr)
    bool _meaningfulGradient(const std::vector<NoisyValue> &grad) const; //check if the gradient is meaningful. i.e. if its values are greater than the statistical errors
    bool _shouldStop(const std::vector<NoisyValue> * grad = nullptr) const; // check for all stopping criteria

    // "Mandatory" logging routines
    // XUpdate should be logged by all optimizers and also GradientInLog if a gradient is used
    void _writeGradientToLog(const std::vector<NoisyValue> &grad) const;
    void _writeXUpdateToLog(const std::vector<double> &xu) const;

    virtual void _findMin() = 0; // child minimization implementation, result to be stored in _last

public:
    explicit NFM(NoisyFunction * targetfun, int max_n_const_values = 20);
    virtual ~NFM() = default;

    // --- Setters
    void setX(int i, double x) { _last.x[i] = x; } // set per element
    void setX(const double x[]); // set via c-style array
    void setX(const std::vector<double> &x) { this->setX(x.data()); } // set via vector
    void setEpsX(double epsx) { _epsx = epsx; }
    void setEpsF(double epsf) { _epsf = epsf; }

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

    double getEpsX() const { return _epsx; }
    double getEpsF() const { return _epsf; }

    // --- Minimization
    void findMin();
};
} // namespace nfm

#endif
