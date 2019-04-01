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
        _epsx(DEFAULT_EPSX), _epsf(DEFAULT_EPSF), _flag_graderr(_flag_gradfun ? _gradfun->hasGradErr() : false),
        _max_n_iterations(0), _max_n_const_values(DEFAULT_MAX_N_CONST), _last(NoisyIOPair(_ndim)) {}

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
        std::vector<double> diff(_last.x.size());
        // deltaX
        std::transform(old.x.begin(), old.x.end(), _last.x.begin(), diff.begin(), std::minus<>()); // diff = old.x-last.x
        _lastDeltaX = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.); // delta = diff.diff
        _lastDeltaX = sqrt(_lastDeltaX); // distance from (original) old.x to last.x
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
    if (_max_n_iterations > 0 && _istep >= _max_n_iterations) {
        LogManager::logString("\nStopping Reason: Maximal iteration count reached.\n");
        return true;
    }
    return false;
}

// --- Protected methods

void NFM::_storeLastValue()
{
    ++_istep; // count step
    this->_writeCurrentXToLog();
    this->_updateDeltas(); // update step deltas

    _old_values.emplace_front(_last);
    while (_old_values.size() > static_cast<size_t>(_max_n_const_values)) {
        _old_values.pop_back();
    }
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


bool NFM::_isGradNoisySmall(const std::vector<NoisyValue> &grad, const bool flag_log) const
{
    if (_flag_graderr && _gradfun->hasGradErr()) {
        for (const NoisyValue &gi : grad) {
            if (gi != 0.) { return false; } // use noisy value overload
        }
        if (flag_log) { LogManager::logString("\nStopping Reason: Gradient is dominated by noise.\n"); }
        return true;
    }
    return false;
}

bool NFM::_shouldStop(const std::vector<NoisyValue> * grad) const
{   // check all stopping criteria
    bool answer = (_isConverged() || !_changedEnough() || _stepLimitReached());
    if (grad != nullptr) {
        answer = answer || _isGradNoisySmall(*grad);
    }
    return answer;
}


// --- Loggers

void NFM::_writeCurrentXToLog() const
{
    if (LogManager::isVerbose()) {
        LogManager::logNoisyIOPair(_last, LogLevel::VERBOSE, "Current position and target value", "x", "f");
    }
    else { LogManager::logNoisyValue(_last.f, LogLevel::NORMAL, "Current target value", "f"); }
}

void NFM::_writeGradientToLog(const std::vector<NoisyValue> &grad) const
{
    LogManager::logNoisyVector(grad, LogLevel::VERBOSE, _flag_graderr, "Raw gradient", "g");
}

void NFM::_writeOldValuesToLog() const
{
    using namespace std;
    stringstream s;
    s << endl << "last values:    ";
    for (const NoisyIOPair &oldv : _old_values) {
        s << oldv.f << "    ";
    }
    s << endl;
    LogManager::logString(s.str());
}

// --- Setters/Getters

void NFM::setX(const double x[])
{
    std::copy(x, x + _ndim, _last.x.data());
}

void NFM::getX(double x[]) const
{
    std::copy(_last.x.data(), _last.x.data() + _ndim, x);
}

// --- findMin

void NFM::findMin()
{
    this->_clearOldValues();
    _istep = 0;
    this->_findMin();
    LogManager::logNoisyIOPair(_last, LogLevel::NORMAL, "Final position and target value");
}
} // namespace nfm