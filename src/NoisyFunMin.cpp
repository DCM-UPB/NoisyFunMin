#include "nfm/NoisyFunMin.hpp"
#include "nfm/LogManager.hpp"

#include <cmath>
#include <sstream>
#include <algorithm>
#include <numeric>


namespace nfm
{

// --- Constructor

NFM::NFM(NoisyFunction * targetfun, const int max_n_const_values):
        _ndim(targetfun->getNDim()), _targetfun(targetfun),
        _gradfun(dynamic_cast<NoisyFunctionWithGradient *>(_targetfun)), _flag_gradfun(_gradfun != nullptr),
        _flag_graderr(_flag_gradfun ? _gradfun->hasGradErr() : false), _epsx(1.e-05), _epsf(0.),
        _last(NoisyIOPair(_ndim)), _max_n_const_values(std::max(1, max_n_const_values)) /*<=1 means check disabled*/
{}

// --- Private methods

bool NFM::_isConverged() const
{
    const auto max_nold = static_cast<size_t>(_max_n_const_values);
    if (max_nold < 2) { return false; }

    if (_old_values.size() >= max_nold) {
        for (auto it = _old_values.begin(); it != _old_values.end(); ++it) {
            if (it != _old_values.begin()) {
                if (!( (*it).f == (*_old_values.begin()).f )) {
                    return false;
                }
            }
        }
        LogManager::logString("\nCost function has stabilised, interrupting minimization procedure.\n");
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
        _lastDeltaF = std::max(0., fabs(_last.f.value - old.f.value) - _last.f.error - old.f.error);
    }
    else { // is first step, initialize deltas
        _lastDeltaX = _epsx; // check will pass
        _lastDeltaF = _epsf; // same
    }
}

bool NFM::_checkDeltas() const
{
    return ((_epsx <= 0. || _lastDeltaX >= _epsx) && (_epsf <= 0. || _lastDeltaF >= _epsf));
}

// --- Protected methods

void NFM::_storeLastValue()
{
    const auto max_nold = static_cast<size_t>(_max_n_const_values);
    this->_writeCurrentXToLog();
    this->_updateDeltas(); // update step deltas

    _old_values.emplace_front(_last);
    if (_old_values.size() > max_nold) {
        _old_values.pop_back();
    }
}

bool NFM::_meaningfulGradient(const std::vector<NoisyValue> &grad) const
{
    if (_flag_graderr) {
        for (auto &gi : grad) {
            if (fabs(gi.value) > gi.error) { return true; }
        }
        LogManager::logString("\nGradient seems to be meaningless, i.e. its error is too large.\n");
        return false;
    }
    return true;
}

bool NFM::_shouldStop(const std::vector<NoisyValue> * grad) const
{   // check all stopping criteria
    bool answer = (_isConverged() || !_checkDeltas());
    if (grad != nullptr) {
        answer = answer || !_meaningfulGradient(*grad);
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

void NFM::_writeXUpdateToLog(const std::vector<double> &xu) const
{
    LogManager::logVector(xu, LogLevel::VERBOSE, "Position update", "u");
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
    s << "equal to first element? ";
    for (auto it = _old_values.begin(); it != _old_values.end(); ++it) {
        if (it != _old_values.begin()) {
            s << ((*it).f == (*_old_values.begin()).f) << "    ";
        }
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
    this->_findMin();
}
} // namespace nfm