#include "nfm/Adam.hpp"

#include "nfm/LogManager.hpp"

#include <algorithm>
#include <cmath>

namespace nfm
{

// --- Constructor

Adam::Adam(NoisyFunctionWithGradient * targetfun, const bool useAveraging,
           const double alpha, const double beta1, const double beta2, const double epsilon):
        NFM(targetfun), _useAveraging(useAveraging), _alpha(alpha), _beta1(beta1), _beta2(beta2), _epsilon(epsilon)
{
    if (!_flag_gradfun) {
        throw std::invalid_argument("[Adam] Adam optimization requires a target function with gradient.");
    }
    // overwrite defaults
    _flag_graderr = false; // don't stop on noisy-low gradients, by default
}

// --- Minimization

void Adam::_findMin()
{
    LogManager::logString("\nBegin Adam::findMin() procedure\n");


    //initialize the gradient & moments
    auto nd = static_cast<size_t>(_ndim);
    std::vector<NoisyValue> grad(nd); // gradient and (unused) error
    std::vector<double> m(nd), v(nd); // moment vectors
    std::vector<double> xavg; // when averaging is enabled, holds the running average
    if (_useAveraging) { xavg.assign(nd, 0.); }

    //begin the minimization loop
    double beta1t = 1.; // stores beta1^t
    double beta2t = 1.; // stores beta2^t
    int iter = 0;
    while (true) {
        ++iter;
        if (LogManager::isLoggingOn()) { // else skip string construction
            LogManager::logString("\nAdam::findMin() Step " + std::to_string(iter) + "\n");
        }

        // compute current gradient and target value
        _last.f = _gradfun->fgrad(_last.x, grad);
        _storeLastValue();
        _writeGradientToLog(grad);
        if (_shouldStop(&grad)) { break; }

        // update factors
        beta1t = beta1t*_beta1; // update beta1 power
        beta2t = beta2t*_beta2; // update beta2 power
        const double afac = _alpha*sqrt(1. - beta2t)/(1. - beta1t);

        // compute the update
        for (int i = 0; i < _ndim; ++i) {
            m[i] = _beta1*m[i] + (1. - _beta1)*grad[i].value; // Update biased first moment
            v[i] = _beta2*v[i] + (1. - _beta2)*grad[i].value*grad[i].value; // Update biased second raw moment

            _last.x[i] -= afac*m[i]/(sqrt(v[i]) + _epsilon); // update _last

            if (_useAveraging) {
                xavg[i] = _beta2*xavg[i] + (1. - _beta2)*_last.x[i];
            }
        }
    }

    if (_useAveraging) { // we need to update _last to the averaged x
        for (int i = 0; i < _ndim; ++i) {
            _last.x[i] = xavg[i]/(1. - beta2t); // bias corrected average
        }
        _last.f = _gradfun->f(_last.x); // evaluate new final function value
    }

    LogManager::logString("\nEnd Adam::findMin() procedure\n");
}
} // namespace nfm
