#include "nfm/DynamicDescent.hpp"

#include "nfm/LogManager.hpp"

#include <numeric>

namespace nfm
{

DynamicDescent::DynamicDescent(NoisyFunctionWithGradient * targetfun, const DDMode ddmode, const bool useAveraging, const double stepSize, const double beta):
        NFM(targetfun), _ddmode(ddmode), _useAveraging(useAveraging), _stepSize(std::max(0., stepSize)), _beta(std::max(0., std::min(1., beta)))
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
    std::vector<double> h; // helper vector used by AdaGrad and RMSProp
    if (_ddmode == DDMode::ADAG || _ddmode == DDMode::RMSP) {
        h.assign(grad.size(), 0.);
    }

    //begin the minimization loop
    int cont = 0;
    while (true) {
        ++cont;
        LogManager::logString("\nDynamicDescent::findMin() Step " + std::to_string(cont) + "\n");

        // compute the gradient and current target
        const bool flag_cont = this->_updateTarget(grad);
        if (!flag_cont) { break; } // we are done

        // find the next position
        this->_findNextX(grad, dx, h);
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

void DynamicDescent::_findNextX(const std::vector<NoisyValue> &grad, std::vector<double> &dx, std::vector<double> &h)
{
    // compute update
    switch (_ddmode) {
    case DDMode::SGDM:
        for (int i = 0; i < _ndim; ++i) {
            dx[i] = _beta*dx[i] - _stepSize*grad[i].value;
        }
        break;
    }
    /*case DDMode::ADAG:
        for
    case DDMode::RMSP:
*/
    //update x
    for (int i = 0; i < _ndim; ++i) {
        _last.x[i] += dx[i];
    }
    // report to the log
    this->_writeXUpdateToLog(dx);
}
} // namespace nfm