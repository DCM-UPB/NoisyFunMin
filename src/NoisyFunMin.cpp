#include "nfm/NoisyFunMin.hpp"

#include "nfm/LogManager.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <sstream>


namespace nfm
{

// --- Constructor

NFM::NFM(NoisyFunction * targetfun):
        _ndim(targetfun->getNDim()), _targetfun(targetfun),
        _gradfun(dynamic_cast<NoisyFunctionWithGradient *>(_targetfun)), _flag_gradfun(_gradfun != nullptr),
        _epsx(DEFAULT_EPSX), _epsf(DEFAULT_EPSF), _flag_gradErrStop(_flag_gradfun ? _gradfun->hasGradErr() : false),
        _max_n_iterations(0), _max_n_const_values(DEFAULT_MAX_N_CONST), _last(_ndim), _grad(_ndim, _gradfun->hasGradErr()) {}

// --- Private methods

bool NFM::_isConverged() const
{
    const auto max_nold = static_cast<size_t>(_max_n_const_values);
    if (max_nold < 2) { return false; } // we need at least two values for this check

    if (_old_values.size() >= max_nold) {
        for (auto &oldp : _old_values) {
            if (oldp.f != _old_values.front().f) {
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
        NoisyIOPair &old = _old_values.front(); // reference to last old value
        // deltaX
        _lastDeltaX = 0.;
        for (int i=0; i<_ndim; ++i) {
            _lastDeltaX += pow(old.x[i]-_last.x[i], 2);
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
    _old_values.emplace_front(_last);
    while (_old_values.size() > static_cast<size_t>(_max_n_const_values)) {
        _old_values.pop_back();
    }

    // call policy
    if (_policy) { _flag_policyStop = _policy(*this, *_targetfun); }

    // count step
    ++_istep;
}

void NFM::_averageOldValues()
{
    std::fill(_last.x.begin(), _last.x.end(), 0.);
    for (auto &oldp : _old_values) {
        std::transform(_last.x.begin(), _last.x.end(), oldp.x.begin(), _last.x.begin(), std::plus<>());
    }
    for (double &x : _last.x) { x /= _old_values.size(); } // get proper averages
    _last.f = _targetfun->f(_last.x); // evaluate final function value
}


bool NFM::_isGradNoisySmall(const bool flag_log) const
{
    if (_flag_gradErrStop && _gradfun->hasGradErr()) {
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
{
    LogManager::logNoisyVector(_grad, LogLevel::VERBOSE, _gradfun->hasGradErr(), "Raw gradient", "g");
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

void NFM::disableStopping()
{ // turn NFM::findMin into an endless loop (unless policy cares for stopping)
    _epsx = 0.;
    _epsf = 0.;
    _flag_gradErrStop = false;
    _max_n_const_values = 1;
    _max_n_iterations = 0;
}

// --- findMin

void NFM::findMin()
{
    this->_clearOldValues();
    _istep = 0;
    _flag_policyStop = false;

    // find minimum
    this->_findMin();
    LogManager::logNoisyIOPair(_last, LogLevel::NORMAL, "Final position and target value");
}
} // namespace nfm
