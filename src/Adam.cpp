#include "nfm/Adam.hpp"

#include "nfm/LogManager.hpp"

#include <algorithm>
#include <cmath>

namespace nfm
{

// --- Minimization

void Adam::findMin()
{
    LogManager::logString("\nBegin Adam::findMin() procedure\n");

    // clear old values
    _clearOldValues();

    //initialize the gradient & moments
    auto nd = static_cast<size_t>(_ndim);
    std::vector<NoisyValue> grad(nd); // gradient and (unused) error
    std::vector<double> m(nd), v(nd); // moment vectors
    std::vector<double> dx(nd); // holds the actual updates for x
    std::vector<double> xavg; // when averaging is enabled, holds the running average
    if (_useAveraging) { xavg.assign(nd, 0.); }

    //begin the minimization loop
    double beta1t = 1.; // stores beta1^t
    double beta2t = 1.; // stores beta2^t
    int step = 0;
    while (true) {
        ++step;

        // compute current gradient and target value
        _last.f = _gradfun->fgrad(_last.x, grad);

        _storeOldValue();
        if (_shouldStop(&grad)) { break; }

        LogManager::logString("\n\nAdam::findMin() Step " + std::to_string(step) + "\n");
        _writeCurrentXToLog();
        _writeGradientToLog(grad);

        beta1t = beta1t*_beta1; // update beta1 power
        beta2t = beta2t*_beta2; // update beat2 power
        const double afac = _alpha*sqrt(1. - beta2t)/(1. - beta1t);

        // compute the update
        for (int i = 0; i < _ndim; ++i) {
            m[i] = _beta1*m[i] + (1. - _beta1)*grad[i].value; // Update biased first moment
            v[i] = _beta2*v[i] + (1. - _beta2)*grad[i].value*grad[i].value; // Update biased second raw moment

            dx[i] = -afac*m[i]/(sqrt(v[i]) + _epsilon); // calculate updates
            _last.x[i] += dx[i]; // update _last

            if (_useAveraging) {
                xavg[i] = _beta2*xavg[i] + (1. - _beta2)*_last.x[i];
            }
        }
        _writeXUpdateToLog(dx);
    }

    if (_useAveraging) { // we need to update _last to the averaged x
        for (int i = 0; i < _ndim; ++i) {
            _last.x[i] = xavg[i]/(1. - beta2t); // bias corrected average
        }
    }

    LogManager::logNoisyIOPair(_last, LogLevel::NORMAL, "Final position and target value");
    LogManager::logString("\nEnd Adam::findMin() procedure\n");
}
} // namespace nfm