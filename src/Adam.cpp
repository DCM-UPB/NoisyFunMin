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

    for (int i=0; i<_ndim; ++i) { // set all to 0
        grad[i] = 0.;
        graderr[i] = 0.;
        m[i] = 0.;
        v[i] = 0.;
        dx[i] = 0.;
    }

    const double afac = _alpha * sqrt(1-_beta2)/(1-_beta1); // helping factor

    //begin the minimization loop
    int step = 0;
    while ( true )
        {
            // compute current gradient and target value
            double newf, newdf;
            this->_gradtargetfun->fgrad(_x->getX(), newf, newdf, grad, graderr);
            _x->setF(newf, newdf);

            if (_shouldStop(grad, graderr)) break;

            log_manager.writeOnLog("\n\nAdam::findMin() Step " + std::to_string(step+1) + "\n");
            _writeCurrentXInLog();
            _writeGradientInLog(grad, graderr);

            // compute the update
            for (int i=0; i<_ndim; ++i) {
                m[i] = _beta1 * m[i] + (1.-_beta1) * grad[i]; // Update biased first moment
                v[i] = _beta2 * v[i] + (1.-_beta2) * grad[i]*grad[i]; // Update biased second raw moment

                dx[i] = - afac * m[i] / (sqrt(v[i]) + _epsilon); // calculate updates
                _x->setX(i, _x->getX(i) + dx[i]); // update _x
            }

            _writeXUpdateInLog(dx);

            step++;
        }

    log_manager.writeNoisyValueInLog(_x, "Final position and target value");
    log_manager.writeOnLog("\nEnd Adam::findMin() procedure\n");
}
