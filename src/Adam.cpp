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
    double grad[_ndim], graderr[_ndim]; // gradient and (unused) error
    double m[_ndim], v[_ndim]; // moment vectors
    double dx[_ndim]; // holds the actual updates for x
    double * xavg = nullptr; // when averaging is enabled, holds the running average
    if (_useAveraging) { xavg = new double[_ndim]; }

    // set all arrays to 0
    std::fill(grad, grad + _ndim, 0.);
    std::fill(graderr, graderr + _ndim, 0.);
    std::fill(m, m + _ndim, 0.);
    std::fill(v, v + _ndim, 0.);
    std::fill(dx, dx + _ndim, 0.);
    if (_useAveraging) { std::fill(xavg, xavg + _ndim, 0.); }

    //begin the minimization loop
    double newf, newdf;
    double beta1t = 1.; // stores beta1^t
    double beta2t = 1.; // stores beta2^t
    int step = 0;
    while (true) {
        ++step;

        // compute current gradient and target value
        this->_gradfun->fgrad(_last->getX(), newf, newdf, grad, graderr);
        _last->setF(newf, newdf);

        _storeOldValue();
        if (_shouldStop(grad, graderr)) { break; }

        LogManager::logString("\n\nAdam::findMin() Step " + std::to_string(step) + "\n");
        _writeCurrentXToLog();
        _writeGradientToLog(grad, graderr);

        beta1t = beta1t*_beta1; // update beta1 power
        beta2t = beta2t*_beta2; // update beat2 power
        const double afac = _alpha*sqrt(1. - beta2t)/(1. - beta1t);

        // compute the update
        for (int i = 0; i < _ndim; ++i) {
            m[i] = _beta1*m[i] + (1. - _beta1)*grad[i]; // Update biased first moment
            v[i] = _beta2*v[i] + (1. - _beta2)*grad[i]*grad[i]; // Update biased second raw moment

            dx[i] = -afac*m[i]/(sqrt(v[i]) + _epsilon); // calculate updates
            _last->setX(i, _last->getX(i) + dx[i]); // update _last

            if (_useAveraging) {
                xavg[i] = _beta2*xavg[i] + (1. - _beta2)*_last->getX(i);
            }
        }
        _writeXUpdateToLog(dx);
    }

    if (_useAveraging) { // we need to update _last to the averaged x
        for (int i = 0; i < _ndim; ++i) {
            _last->setX(i, xavg[i]/(1. - beta2t)); // bias corrected average
        }
        this->_gradfun->fgrad(_last->getX(), newf, newdf, grad, graderr);
        _last->setF(newf, newdf);

        delete[] xavg;
    }

    LogManager::logNoisyIOPair(_last, 1, "Final position and target value");
    LogManager::logString("\nEnd Adam::findMin() procedure\n");
}
} // namespace nfm