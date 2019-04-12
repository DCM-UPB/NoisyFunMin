#include "nfm/NoisyFunMin.hpp"

#include "nfm/LogManager.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>


namespace nfm
{

// --- Constructor

NFM::NFM(const int ndim, const bool needsGrad):
        _ndim(ndim), _flag_needsGrad(needsGrad), _last(_ndim), _grad(_ndim), _flag_gradErrStop(needsGrad /*default*/)
{
    _old_values.reserve(static_cast<size_t>(_max_n_const_values));
}


// --- Private methods

bool NFM::_isConverged() const
{
    const auto max_nold = static_cast<size_t>(_max_n_const_values);
    if (max_nold < 2) { return false; } // we need at least two values for this check

    if (_old_values.size() >= max_nold) {
        for (const auto &oldp : _old_values.vec()) {
            if (oldp.f != _old_values.back().f) {
                return false;
            }
        }
        LogManager::logString("\nStopping Reason: Target function list has stabilized.\n");
        return true;
    }

    return false;
}

void NFM::_updateDeltas()
{
    if (!_old_values.empty()) {
        const NoisyIOPair &old = _old_values.back(); // reference to last old value
        // deltaX
        _lastDeltaX = 0.;
        for (int i = 0; i < _ndim; ++i) {
            _lastDeltaX += pow(old.x[i] - _last.x[i], 2);
        }
        _lastDeltaX = sqrt(_lastDeltaX);
        // deltaF
        _lastDeltaF = std::max(0., _last.f.minDist(old.f));
    }
    else { // is first step, initialize deltas
        _lastDeltaX = _epsx; // check will pass
        _lastDeltaF = _epsf; // same
    }
}

bool NFM::_changedEnough() const
{
    if (_epsx > 0. && _lastDeltaX < _epsx) {
        LogManager::logString("\nStopping Reason: Position did not change enough.\n");
        return false;
    }
    if (_epsf > 0. && _lastDeltaF < _epsf) {
        LogManager::logString("\nStopping Reason: Target function did not change enough.\n");
        return false;
    }
    return true;
}

bool NFM::_stepLimitReached() const
{
    if (_max_n_iterations > 0 && _istep > _max_n_iterations) {
        LogManager::logString("\nStopping Reason: Maximal iteration count reached.\n");
        return true;
    }
    return false;
}

// --- Protected methods

void NFM::_storeLastValue()
{
    this->_writeCurrentXToLog();
    this->_updateDeltas(); // changes in x and f

    // update old value list
    _old_values.push_back(_last); // oldest element will be deleted (if full)

    // call policy
    if (_policy) { _flag_policyStop = _policy(*this, *_targetfun); }

    // count step
    ++_istep;
}

void NFM::_averageOldValues()
{
    std::fill(_last.x.begin(), _last.x.end(), 0.);
    for (const auto &oldp : _old_values.vec()) {
        std::transform(_last.x.begin(), _last.x.end(), oldp.x.begin(), _last.x.begin(), std::plus<>());
    }
    for (double &x : _last.x) { x /= _old_values.size(); } // get proper averages
    _last.f = _targetfun->f(_last.x); // evaluate final function value
}


bool NFM::_isGradNoisySmall(const bool flag_log) const
{
    if (_flag_gradErrStop && this->hasGradErr()) {
        if (_grad > 0.) { return false; } // use overload
        if (flag_log) { LogManager::logString("\nStopping Reason: Gradient is dominated by noise.\n"); }
        return true;
    }
    return false;
}

bool NFM::_shouldStop() const
{   // check all stopping criteria
    if (_flag_policyStop) {
        LogManager::logString("\nStopping Reason: User provided policy.\n");
        return true;
    }
    return (_isConverged() || !_changedEnough() || _stepLimitReached() || _isGradNoisySmall());
}


// --- Loggers

void NFM::_writeCurrentXToLog() const
{
    if (LogManager::isVerbose()) {
        LogManager::logNoisyIOPair(_last, LogLevel::VERBOSE, "Current position and target value", "x", "f");
    }
    else { LogManager::logNoisyValue(_last.f, LogLevel::NORMAL, "Current target value", "f"); }
}

void NFM::_writeGradientToLog() const
{   // !! Derived: Use only if NFM constructor is called with needsErr = true
    LogManager::logNoisyVector(_grad, LogLevel::VERBOSE, this->hasGradErr(), "Raw gradient", "g");
}

// --- Setters/Getters

void NFM::setX(const double x[])
{
    std::copy(x, x + _ndim, _last.x.data());
}

void NFM::setX(const std::vector<double> &x)
{
    if (x.size() == _last.x.size()) {
        _last.x = x;
    }
    else {
        throw std::invalid_argument("[NFM::setX] Passed vector length didn't match NFM's number of dimensions.");
    }
}

void NFM::getX(double x[]) const
{
    std::copy(_last.x.data(), _last.x.data() + _ndim, x);
}

void NFM::setMaxNConstValues(int maxn_const_values)
{
    _max_n_const_values = std::max(1, maxn_const_values);
    _old_values.set_cap(static_cast<size_t>(_max_n_const_values));
}

void NFM::disableStopping()
{ // turn NFM::findMin into an endless loop (unless policy cares for stopping)
    _epsx = 0.;
    _epsf = 0.;
    _flag_gradErrStop = false;
    this->setMaxNConstValues(1);
    _max_n_iterations = 0;
}

// --- findMin

NoisyIOPair NFM::findMin(NoisyFunction &targetfun)
{
    if (targetfun.getNDim() != this->getNDim()) {
        throw std::invalid_argument("[NFM] Passed target function's number of inputs is not equal to NFM's number of dimensions.");
    }

    // setup target function
    _targetfun = &targetfun; // we keep a pointer during findMin()
    _gradfun = dynamic_cast<NoisyFunctionWithGradient *>(_targetfun); // we do this single dynamic cast to check for gradient functions

    if (this->needsGrad()) { // setup gradient case
        if (_gradfun != nullptr) {
            _flag_validGrad = true; // gradient values will be calulcated and stored
            _flag_validGradErr = _gradfun->hasGradErr(); // for the errors it also depends on the function's settings
        }
        else {
            throw std::invalid_argument("[NFM] The optimizer requires gradients, but the target function doesn't provide them.");
        }
    }
    else { // gradients will not be calculated
        _flag_validGrad = false;
        _flag_validGradErr = false;
    }

    // for consistency reset some values
    _last.f.zero();
    _grad.zero();
    _old_values.clear();
    _lastDeltaX = 0.;
    _lastDeltaF = 0.;
    _istep = 0;
    _flag_policyStop = false;

    // find minimum
    this->_findMin();
    LogManager::logNoisyIOPair(_last, LogLevel::NORMAL, "Final position and target value");

    // reset temporary pointers
    _targetfun = nullptr;
    _gradfun = nullptr;

    return _last;
}

NoisyIOPair NFM::findMin(NoisyFunction &targetFun, const std::vector<double> &x0)
{
    this->setX(x0);
    return this->findMin(targetFun);
}

NoisyIOPair NFM::findMin(NoisyFunction &targetFun, const double x0[])
{
    this->setX(x0);
    return this->findMin(targetFun);
}
} // namespace nfm
