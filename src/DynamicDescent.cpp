#include "nfm/DynamicDescent.hpp"

#include "nfm/LogManager.hpp"

#include <numeric>
#include <cmath>

namespace nfm
{

DynamicDescent::DynamicDescent(NoisyFunctionWithGradient * targetfun, const DDMode ddmode, const bool useAveraging, const double stepSize):
        NFM(targetfun), _ddmode(ddmode), _useAveraging(useAveraging), _stepSize(std::max(0., stepSize))
{
    if (!_flag_gradfun) {
        throw std::invalid_argument("[DynamicDescent] Dynamic Descent optimization requires a target function with gradient.");
    }
    // overwrite defaults
    this->setGradErrStop(false); // don't stop on noisy-low gradients, by default
}

// --- Minimization

void DynamicDescent::_findMin()
{
    LogManager::logString("\nBegin DynamicDescent::findMin() procedure\n");

    // vectors used for SGD updates
    std::vector<double> v(_grad.size()); // helper vector used by all methods
    std::vector<double> w; // only used by AdaDelta
    if (_ddmode == DDMode::ADAD) {
        w = v; // v is already all 0
    }

    //begin the minimization loop
    int iter = 0;
    while (true) {
        ++iter;
        if (LogManager::isLoggingOn()) { // else skip string construction
            LogManager::logString("\nDynamicDescent::findMin() Step " + std::to_string(iter) + "\n");
        }

        // compute the gradient and current target
        this->_updateTarget();
        if (this->_shouldStop()) { break; } // we are done

        // find the next position
        this->_findNextX(iter, v, w);
    }

    if (_useAveraging) { // calculate the old value average as end result
        this->_averageOldValues(); // perform average and store it in last
    }

    LogManager::logString("\nEnd DynamicDescent::findMin() procedure\n");
}

// --- Internal methods

void DynamicDescent::_updateTarget()
{
    _last.f = _gradfun->fgrad(_last.x, _grad);
    this->_storeLastValue();
    this->_writeGradientToLog();
}

void DynamicDescent::_findNextX(const int iter, std::vector<double> &v, std::vector<double> &w)
{
    const std::vector<double> &gradv = _grad.val; // we need only the values

    // compute update
    switch (_ddmode) {
    case DDMode::SGDM:
        for (int i = 0; i < _ndim; ++i) {
            v[i] = _beta*v[i] + _stepSize*gradv[i];
            _last.x[i] += v[i];
        }
        break;

    case DDMode::ADAG:
        for (int i = 0; i < _ndim; ++i) {
            v[i] += gradv[i]*gradv[i];
            _last.x[i] += _stepSize/(sqrt(v[i]) + _epsilon)*gradv[i];
        }
        break;

    case DDMode::ADAD:
        if (iter > 1) {
            for (int i = 0; i < _ndim; ++i) {
                v[i] = _beta*v[i] + (1. - _beta)*(gradv[i]*gradv[i]);
                const double dx = gradv[i]*(sqrt(w[i]) + _epsilon)/(sqrt(v[i]) + _epsilon);
                _last.x[i] += dx;
                w[i] = _beta*w[i] + (1. - _beta)*(dx*dx);
            }
        }
        else { // first step
            for (int i = 0; i < _ndim; ++i) {
                v[i] = gradv[i]*gradv[i];
                const double dx = _stepSize*gradv[i]; // initially we use the stepSize
                _last.x[i] += dx;
                w[i] = dx*dx;
            }
        }
        break;

    case DDMode::RMSP: // standard RMSProp with first step as grad desc
        if (iter > 1) {
            for (int i = 0; i < _ndim; ++i) {
                v[i] = _beta*v[i] + (1. - _beta)*(gradv[i]*gradv[i]);
                _last.x[i] += _stepSize*gradv[i]/(sqrt(v[i]) + _epsilon);
            }
        }
        else { // first step
            for (int i = 0; i < _ndim; ++i) {
                v[i] = gradv[i]*gradv[i];
                _last.x[i] += _stepSize*gradv[i];
            }
        }
        break;

    case DDMode::NEST: // Bengio update with first step as grad desc
        if (iter > 1) {
            for (int i = 0; i < _ndim; ++i) {
                _last.x[i] += _beta*_beta*v[i] + (1. + _beta)*_stepSize*gradv[i];
                v[i] = _beta*v[i] + _stepSize*gradv[i];
            }
        }
        else { // first step
            for (int i = 0; i < _ndim; ++i) {
                _last.x[i] += _stepSize*gradv[i];
            }
        }
    }
}
} // namespace nfm
