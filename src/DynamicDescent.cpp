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
    std::vector<double> v(grad.size()); // helper vector used by all methods
    std::vector<double> w; // only used by AdaDelta
    if (_ddmode == DDMode::ADAD) {
        w = v; // v is already all 0
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
        this->_findNextX(iter, grad, v, w);
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

void DynamicDescent::_findNextX(const int iter, const std::vector<NoisyValue> &grad, std::vector<double> &v, std::vector<double> &w)
{
    // compute update
    switch (_ddmode) {
    case DDMode::SGDM:
        for (int i = 0; i < _ndim; ++i) {
            v[i] = _beta*v[i] - _stepSize*grad[i].value;
            _last.x[i] += v[i];
        }
        break;

    case DDMode::ADAG:
        for (int i = 0; i < _ndim; ++i) {
            v[i] += grad[i].value*grad[i].value;
            _last.x[i] -= _stepSize/(sqrt(v[i]) + _epsilon)*grad[i].value;
        }
        break;

    case DDMode::ADAD:
        if (iter > 1) {
            for (int i = 0; i < _ndim; ++i) {
                v[i] = _beta*v[i] + (1. - _beta)*(grad[i].value*grad[i].value);
                const double dx = -grad[i].value*(sqrt(w[i]) + _epsilon)/(sqrt(v[i]) + _epsilon);
                _last.x[i] += dx;
                w[i] = _beta*w[i] + (1. - _beta)*(dx*dx);
            }
        }
        else { // first step
            for (int i = 0; i < _ndim; ++i) {
                v[i] = grad[i].value*grad[i].value;
                const double dx = -_stepSize*grad[i].value; // initially we use the stepSize
                _last.x[i] += dx;
                w[i] = dx*dx;
            }
        }
        break;

    case DDMode::RMSP: // standard RMSProp with first step as grad desc
        if (iter > 1) {
            for (int i = 0; i < _ndim; ++i) {
                v[i] = _beta*v[i] + (1. - _beta)*(grad[i].value*grad[i].value);
                _last.x[i] -= _stepSize*grad[i].value/(sqrt(v[i]) + _epsilon);
            }
        }
        else { // first step
            for (int i = 0; i < _ndim; ++i) {
                v[i] = grad[i].value*grad[i].value;
                _last.x[i] -= _stepSize*grad[i].value;
            }
        }
        break;

    case DDMode::NEST: // Bengio update with first step as grad desc
        if (iter > 1) {
            for (int i = 0; i < _ndim; ++i) {
                _last.x[i] += _beta*_beta*v[i] - (1. + _beta)*_stepSize*grad[i].value;
                v[i] = _beta*v[i] - _stepSize*grad[i].value;
            }
        }
        else { // first step
            for (int i = 0; i < _ndim; ++i) {
                _last.x[i] -= _stepSize*grad[i].value;
            }
        }
    }
}
} // namespace nfm
