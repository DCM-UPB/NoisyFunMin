#ifndef NFM_NOISYFUNMIN_HPP
#define NFM_NOISYFUNMIN_HPP

#include <functional>

#include "nfm/NoisyFunction.hpp"
#include "nfm/PushBackBuffer.hpp"

namespace nfm
{

class NFM
{
protected:
    // Consts
    const int _ndim;  //dimensionality of the space where the target function is embedded
    const bool _flag_needsGrad; // does the optimization method use gradients?

    // Temporary pointers valid during findMin()
    NoisyFunction * _targetfun{};  // target function to minimize
    NoisyFunctionWithGradient * _gradfun{};  // will be null if _targetfun is not NoisyFunctionWithGradient

    // Member to be updated by child
    NoisyIOPair _last; // last position and its function value (will be 0 after construction)
    NoisyGradient _grad; // last noisy gradient (will be all 0 if current/last run didn't use gradients)

private: // set/called directly by base class only
    PushBackBuffer<NoisyIOPair> _old_values; // list of previous target values and positions

    // Stopping Criteria (childs should use get/set methods, but may change defaults)
    bool _flag_gradErrStop; // should we consider gradient errors for stopping? (if targetfun supports it)
    double _epsx = 1.e-5; // changes in the position x smaller than this value will stop the minimization (if 0, disabled)
    double _epsf = 0.; // changes in the function smaller than this value will stop the minimization (if 0, disabled)
    int _max_n_iterations = 0; // hard stop after this amount of iterations (if 0, disabled (the default!))
    int _max_n_const_values = 20; // stop after this number of target values have been constant within error bounds (if <= 1, disabled)

    // other
    double _lastDeltaX{}; // change in x by last step (updated on storeLastValue)
    double _lastDeltaF{}; // change in f by last step (updated on storeLastalue)
    int _istep{}; // counts the calls to _storeLastValue()
    bool _flag_policyStop{}; // did the user policy dictate stopping?
    bool _flag_validGrad{}; // are the values in _grad meaningful (i.e. not default 0)
    bool _flag_validGradErr{}; // if the targetfun doesn't provide gradient errors, this stays false

    // Optional user provided policy function, called after every iteration.
    // May manipulate the NFM and target function (passed as their base types).
    // If necessary, upcast them to their known true type. Must return true
    // for NFM to continue, else NFM will stop at the next shouldStop() check.
    std::function<bool(NFM &, NoisyFunction &)> _policy{};

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
    bool _isGradNoisySmall(bool flag_log = true) const; //check if any gradient element is greater than its statistical error
    bool _shouldStop() const; // check for all stopping criteria

    // "Mandatory" logging routines
    // If a gradient is used, it should be logged after it is calculated
    void _writeGradientToLog() const;

    // TO BE IMPLEMENTED
    virtual void _findMin() = 0; // minimization implementation (called in findMin(), result to be stored in _last

public:
    NFM(int ndim, bool needsGrad);
    virtual ~NFM() = default;

    // --- Setters

    // You may pre-set the initial position (or pass x0 on findMin())
    void setX(int i, double x) { _last.x[i] = x; } // set per element
    void setX(const double x[]); // set via c-style array
    void setX(const std::vector<double> &x); // set via vector (must be ndim length)

    // stopping conditions
    void setEpsX(double epsx) { _epsx = epsx; }
    void setEpsF(double epsf) { _epsf = epsf; }
    void setGradErrStop(bool flag_gradErrStop) { _flag_gradErrStop = flag_gradErrStop; }
    void setMaxNIterations(int maxn_iterations) { _max_n_iterations = std::max(0, maxn_iterations); }
    void setMaxNConstValues(int maxn_const_values);
    void disableStopping(); // WILL TURN OFF ALL STOPPING CRITERIA (except user policy)

    // Set an own policy function which may manipulate NFM and target function on each step.
    // It will always get called after a new position pair has been stored.
    void setPolicy(const std::function<bool(NFM &, NoisyFunction &)> &policy) { _policy = policy; }
    void clearPolicy() { _policy = nullptr; } // set empty policy

    // --- Getters

    int getNDim() const { return _ndim; }

    // Last Position
    double getX(int i) const { return _last.x[i]; }; // elementary get
    void getX(double x[]) const; // get via c-style array
    void getX(std::vector<double> &x) const { this->getX(x.data()); } // get via passed vector
    const std::vector<double> &getX() const { return _last.x; } // get const ref

    // Last Value or Value/Position Pair
    double getF() const { return _last.f.val; }
    double getDf() const { return _last.f.err; }
    NoisyValue getFDf() const { return _last.f; }
    const NoisyIOPair &getLast() const { return _last; }

    // Gradient (will be all 0 if hasGrad() == false)
    const NoisyGradient &getGrad() const { return _grad; }
    bool hasGrad() const { return _flag_validGrad; } // is getGrad().val meaningful?
    bool hasGradErr() const { return _flag_validGradErr; } // is getGrad().err meaningful?
    bool needsGrad() const { return _flag_needsGrad; } // does the derived optimizer require gradients?

    // Other last values
    const PushBackBuffer<NoisyIOPair> &getOldValues() const { return _old_values; }
    double getDeltaX() const { return _lastDeltaX; }
    double getDeltaF() const { return _lastDeltaF; }
    double getIter() const { return _istep; }

    // Stopping
    double getEpsX() const { return _epsx; }
    double getEpsF() const { return _epsf; }
    bool getGradErrStop() const { return _flag_gradErrStop; }
    int getMaxNIterations() const { return _max_n_iterations; }
    int getMaxNConstValues() const { return _max_n_const_values; }


    // When in your use case (for whatever reason) it can happen that you access
    // a NFM object while it is running the findMin() method, this may be used to check.
    bool isRunning() const { return _targetfun != nullptr; }


    // --- Minimization method

    // Minimize a target function (must be a NoisyFunctionWithGradient when needsGrad()).
    // The initial position will be the last internal X position (0 after construction).
    // Returns the optimal position and function value in a NoisyIOPair.
    NoisyIOPair findMin(NoisyFunction &targetFun); // will throw if wrong ndim or gradient requirement not matched

    // Optionally provide different initial positions x0
    NoisyIOPair findMin(NoisyFunction &targetFun, const std::vector<double> &x0); // will throw on wrong size
    NoisyIOPair findMin(NoisyFunction &targetFun, const double x0[]); // c array version (no size check)
};
} // namespace nfm

#endif
