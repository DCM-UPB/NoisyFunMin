#include "Adam.hpp"

#include "LogNFM.hpp"

#include <cmath>


// --- Minimization

void Adam::findMin(){
    NFMLogManager log_manager = NFMLogManager();
    log_manager.writeOnLog("\nBegin Adam::findMin() procedure\n");

    // clear old values
    _clearOldValues();

    //initialize the gradient & moments
    double grad[_ndim], graderr[_ndim]; // gradient and (unused) error
    double m[_ndim], v[_ndim]; // moment vectors
    double dx[_ndim]; // holds the actual updates for x
    double * xavg; // when averaging is enabled, holds the running average
    if (_useAveraging) xavg = new double[_ndim];

    for (int i=0; i<_ndim; ++i) { // set all to 0
        grad[i] = 0.;
        graderr[i] = 0.;
        m[i] = 0.;
        v[i] = 0.;
        dx[i] = 0.;
        if (_useAveraging) xavg[i] = 0.;
    }

    //begin the minimization loop
    double newf, newdf;
    double beta1t = 1.; // stores beta1^t
    double beta2t = 1.; // stores beta2^t
    int step = 0;
    while ( true )
        {
            ++step;

            // compute current gradient and target value
            this->_gradtargetfun->fgrad(_x->getX(), newf, newdf, grad, graderr);
            _x->setF(newf, newdf);

            if (_shouldStop(grad, graderr)) break;

            log_manager.writeOnLog("\n\nAdam::findMin() Step " + std::to_string(step) + "\n");
            _writeCurrentXInLog();
            _writeGradientInLog(grad, graderr);

            // compute the update
            for (int i=0; i<_ndim; ++i) {
                m[i] = _beta1 * m[i] + (1.-_beta1) * grad[i]; // Update biased first moment
                v[i] = _beta2 * v[i] + (1.-_beta2) * grad[i]*grad[i]; // Update biased second raw moment

                beta1t = beta1t * _beta1; // update beta1 power
                beta2t = beta2t * _beta2; // update beat2 power

                dx[i] = - _alpha * sqrt(1.-beta2t) / (1.-beta1t) * m[i] / (sqrt(v[i]) + _epsilon); // calculate updates
                _x->setX(i, _x->getX(i) + dx[i]); // update _x

                if (_useAveraging) {
                    xavg[i] = _beta2 * xavg[i] + (1.-_beta2) * _x->getX(i);
                    xavg[i] /= (1.-beta2t);
                }
            }
            _writeXUpdateInLog(dx);
        }

    if (_useAveraging) { // we need to update _x to the averaged x
        for (int i=0; i<_ndim; ++i) _x->setX(i, xavg[i]);
        this->_gradtargetfun->fgrad(_x->getX(), newf, newdf, grad, graderr);
        _x->setF(newf, newdf);
    }

    log_manager.writeNoisyValueInLog(_x, "Final position and target value");
    log_manager.writeOnLog("\nEnd Adam::findMin() procedure\n");
}
