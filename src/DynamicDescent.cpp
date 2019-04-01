#include "nfm/DynamicDescent.hpp"

#include "nfm/LogManager.hpp"

#include <numeric>

namespace nfm
{

DynamicDescent::DynamicDescent(NoisyFunctionWithGradient * targetfun, const DDMode ddmode, const bool useAveraging, const double stepSize, const double beta, const double epsilon):
        NFM(targetfun), _ddmode(ddmode), _useAveraging(useAveraging), _stepSize(std::max(0., stepSize)), _beta(std::max(0., std::min(1., beta))), _epsilon(std::max(0., epsilon))
{
    if (!_flag_gradfun) {
        throw std::invalid_argument("[DynamicDescent] Dynamic Descent optimization requires a target function with gradient.");
    }
    // overwrite defaults
    _flag_graderr = false; // don't stop on noisy-low gradients, by default
}

// --- Minimization

void DynamicDescent::_findMin()
{
    LogManager::logString("\nBegin DynamicDescent::findMin() procedure\n");

    // vector to hold the gradient and (possibly) error
    std::vector<NoisyValue> grad(static_cast<size_t>(_ndim));

    // holds previous old update (init to 0)
    std::vector<double> dx(grad.size()); // stores momentum updates
    std::vector<double> h; // helper vector used by all but SGDM
    if (_ddmode != DDMode::SGDM) {
        h = dx; // dx is all 0
    }

    //begin the minimization loop
    int iter = 0;
    while (true) {
        ++iter;
        LogManager::logString("\nDynamicDescent::findMin() Step " + std::to_string(iter) + "\n");

        // compute the gradient and current target
        const bool flag_cont = this->_updateTarget(grad);
        if (!flag_cont) { break; } // we are done

        // find the next position
        this->_findNextX(grad, dx, h, iter);
    }

    if (_useAveraging) { // calculate the old value average as end result
        this->_averageOldValues(); // perform average and store it in last
    }
    LogManager::logString("\nEnd DynamicDescent::findMin() procedure\n");
}

// --- Internal methods

bool DynamicDescent::_updateTarget(std::vector<NoisyValue> &grad)
{
    _last.f = _gradfun->fgrad(_last.x, grad);
    this->_storeLastValue();
    this->_writeGradientToLog(grad);
    return !this->_shouldStop(&grad);
}

void DynamicDescent::_findNextX(const std::vector<NoisyValue> &grad, std::vector<double> &dx, std::vector<double> &h, const int iter)
{
    // compute update
    switch (_ddmode) {
    case DDMode::SGDM:
        for (int i = 0; i < _ndim; ++i) {
            dx[i] = _beta*dx[i] - _stepSize*grad[i].value;
        }
        break;
    case DDMode::ADAG:
        for (int i = 0; i < _ndim; ++i) {
            h[i] += grad[i].value*grad[i].value;
            dx[i] = -_stepSize/(sqrt(h[i]) + _epsilon)*grad[i].value;
        }
        break;
    case DDMode::RMSP: // standard RMSProp with first step as grad desc
        for (int i = 0; i < _ndim; ++i) {
            if (iter > 1) {
                h[i] = _beta*h[i] + (1. - _beta)*(grad[i].value*grad[i].value);
                dx[i] = -_stepSize*grad[i].value/(sqrt(h[i]) + _epsilon);
            }
            else { // first step
                h[i] = grad[i].value*grad[i].value;
                dx[i] = -_stepSize*grad[i].value;
            }
        }
        break;
    case DDMode::NEST: // Bengio-like update with first step as grad desc
        for (int i = 0; i < _ndim; ++i) {
            if (iter > 1) {
                dx[i] = _beta*_beta*h[i] - (1. + _beta)*_stepSize*grad[i].value;
                h[i] = _beta*h[i] - _stepSize*grad[i].value;
            }
            else { // first step
                dx[i] = -_stepSize*grad[i].value;
            }
        }
    }

    //update x
    for (int i = 0; i < _ndim; ++i) {
        _last.x[i] += dx[i];
    }
    // report to the log
    this->_writeXUpdateToLog(dx);
}
} // namespace nfm