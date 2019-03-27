#include "nfm/NoisyFunMin.hpp"
#include "nfm/LogManager.hpp"

#include <cmath>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <nfm/NoisyFunMin.hpp>


namespace nfm
{

// --- Protected methods


void NFM::_storeOldValue()
{
    const auto max_nold = static_cast<size_t>(_max_n_const_values);
    if (max_nold  > 0) {
        _old_values.emplace_front(NoisyValue(_last.f));

        if (_old_values.size() > max_nold) {
            _old_values.pop_back();
        }
    }
}


bool NFM::_isConverged() const
{
    const auto max_nold = static_cast<size_t>(_max_n_const_values);

    if (max_nold < 1) { return false; }

    if (_old_values.size() == max_nold) {
        for (auto it = _old_values.begin(); it != _old_values.end(); ++it) {
            if (it != _old_values.begin()) {
                if (!(*it == *_old_values.begin())) {
                    return false;
                }
            }
        }
        LogManager::logString("\nCost function has stabilised, interrupting minimization procedure.\n");
        return true;
    }

    return false;
}


bool NFM::_meaningfulGradient(const std::vector<NoisyValue> &grad) const
{
    if (_flag_graderr) {
        for (auto &gi : grad) {
            if (fabs(gi.value) > gi.error) { return true; }
        }
    }
    else {
        return true;
    }

    LogManager::logString("\nGradient seems to be meaningless, i.e. its error is too large.\n");
    return false;
}

bool NFM::_shouldStop(const std::vector<NoisyValue> * grad) const
{
    if (grad != nullptr) {
        return _isConverged() || !_meaningfulGradient(*grad);
    }
    else {
        return _isConverged();
    }
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

void NFM::_writeOldValuesToLog()
{
    using namespace std;

    stringstream s;
    s << endl << "last values:    ";
    for (NoisyValue old_value : _old_values) {
        s << old_value << "    ";
    }
    s << endl;
    s << "equal to first element? ";
    for (auto it = _old_values.begin(); it != _old_values.end(); ++it) {
        if (it != _old_values.begin()) {
            s << ((*it) == (*_old_values.begin())) << "    ";
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

// --- Constructor and destructor

NFM::NFM(NoisyFunction * targetfun, const int max_n_const_values):
        _ndim(targetfun->getNDim()), _targetfun(targetfun),
        _gradfun(dynamic_cast<NoisyFunctionWithGradient *>(_targetfun)), _flag_gradfun(_gradfun != nullptr),
        _flag_graderr(_flag_gradfun ? _gradfun->hasGradErr() : false), _max_n_const_values(max_n_const_values),
        _last(NoisyIOPair(_ndim)), _epsf(0.), _epsx(0.)
        {}
} // namespace nfm