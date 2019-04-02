#ifndef NFM_NOISYFUNMIN_HPP
#define NFM_NOISYFUNMIN_HPP

#include <list>
#include <functional>

#include "nfm/NoisyFunction.hpp"
#include "nfm/NoisyValue.hpp"

namespace nfm
{

class NFM
{
protected:
    // Default values for general NFM parameters (childs may use different defaults)
    static constexpr int DEFAULT_MAX_N_CONST = 20; // stop on maximal number of noisy-constant values
    static constexpr double DEFAULT_EPSX = 1.e-5; // stop on minimal position change
    static constexpr double DEFAULT_EPSF = 0.; // stop on minimal target function change

    // Const members
    const int _ndim;  //dimensionality of the space where the target function is embedded
    NoisyFunction * const _targetfun;  //target function to minimize
    NoisyFunctionWithGradient * const _gradfun;  //gradient of the target function
    const bool _flag_gradfun; // has the gradient been provided?

    // Stopping Tolerances (may also be used by child to auto-set own tolerances, e.g. for line search)
    double _epsx; // changes in the position x smaller than this value will stop the minimization (if 0, disabled)
    double _epsf; // changes in the function smaller than this value will stop the minimization (if 0, disabled)
    bool _flag_gradErrStop; // should we consider gradient errors for stopping? (if targetfun supports it)
    int _max_n_iterations; // hard stop after this amount of iterations (if 0, disabled (the default!))
    int _max_n_const_values; // stop after this number of target values have been constant within error bounds (if <= 1, disabled)

    // Member to be updated by child
    NoisyIOPair _last; // last position and its function value

private: // set/called by base class only
    std::list<NoisyIOPair> _old_values{}; // list of previous target values and positions
    double _lastDeltaX{}; // change in x by last step (updated on storeLastValue)
    double _lastDeltaF{}; // change in f by last step (updated on storeLastalue)
    int _istep{}; // counts the calls to _storeLastValue()

    std::function< void(NFM &, NoisyFunction &) > _policy{}; // optional user provided policy function

    void _clearOldValues() { _old_values.clear(); } // reset old values list
    bool _isConverged() const; // check if the target function has stabilized
    void _updateDeltas(); // calculate deltaX and deltaF between _last and _old_values.front()
    bool _changedEnough() const; // check deltas against epsx and epsf
    bool _stepLimitReached() const; // is the set maximum amount of iteration reached

    void _writeCurrentXToLog() const; // write current x on log on storeLastValue

protected: // Protected methods for child optimizers
    // use this after every position&function update
    void _storeLastValue(); // store last value in old values list (updates deltax/deltaf)

    // use this before terminating, if averaging is desired
    void _averageOldValues(); // compute average x of old value list, store it with the corresponding function value in last

    // check stopping criteria (shouldStop contains meaningfulGradient, if grad!=nullptr)
    bool _isGradNoisySmall(const std::vector<NoisyValue> &grad, bool flag_log = true) const; //check if any gradient element is greater than its statistical error
    bool _shouldStop(const std::vector<NoisyValue> * grad = nullptr) const; // check for all stopping criteria

    // "Mandatory" logging routines
    // If a gradient is used, it should be logged at every step
    void _writeGradientToLog(const std::vector<NoisyValue> &grad) const;

    // TO BE IMPLEMENTED
    virtual void _findMin() = 0; // minimization implementation, result to be stored in _last

public:
    explicit NFM(NoisyFunction * targetfun);
    virtual ~NFM() = default;

    // --- Setters

    void setX(int i, double x) { _last.x[i] = x; } // set per element
    void setX(const double x[]); // set via c-style array
    void setX(const std::vector<double> &x); // set via vector (must be ndim length)

    void setEpsX(double epsx) { _epsx = epsx; }
    void setEpsF(double epsf) { _epsf = epsf; }
    void setGradErrStop(bool flag_gradErrStop) { _flag_gradErrStop = flag_gradErrStop; }
    void setMaxNIterations(int maxn_iterations) { _max_n_iterations = std::max(0, maxn_iterations); }
    void setMaxNConstValues(int maxn_const_values) { _max_n_const_values = std::max(1, maxn_const_values); }

    // Set an own policy function which may manipulate NFM and target function on each step.
    // It will always get called after a new position pair has been stored.
    void setPolicy(const std::function<void(NFM &, NoisyFunction &)> &policy) { _policy = policy; }
    void clearPolicy() { _policy = nullptr; } // set empty policy

    // --- Getters

    int getNDim() const { return _ndim; }
    NoisyFunction * getTargetFun() const { return _targetfun; }
    NoisyFunctionWithGradient * getGradientFun() const { return _gradfun; }

    double getX(int i) const { return _last.x[i]; }; // elementary get
    void getX(double x[]) const; // get via c-style array
    void getX(std::vector<double> &x) const { this->getX(x.data()); } // get via passed vector
    const std::vector<double> &getX() const { return _last.x; } // get const ref
    double getF() const { return _last.f.value; }
    double getDf() const { return _last.f.error; }
    NoisyValue getFDf() const { return _last.f; }
    const NoisyIOPair &getLast() const { return _last; }

    double getEpsX() const { return _epsx; }
    double getEpsF() const { return _epsf; }
    bool getGradErrStop() const { return _flag_gradErrStop; }
    int getMaxNIterations() const { return _max_n_iterations; }
    int getMaxNConstValues() const { return _max_n_const_values; }

    const std::list<NoisyIOPair> & getOldValues() const { return _old_values; }
    double getDeltaX() const { return _lastDeltaX; }
    double getDeltaF() const { return _lastDeltaF; }
    double getIter() const { return _istep; }


    // --- Minimization
    void findMin();
};
} // namespace nfm

#endif
