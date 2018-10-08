#include "Adam.hpp"

#include "LogNFM.hpp"

#include <cmath>


// --- Minimization

void Adam::findMin(){
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

    NFMLogManager log_manager = NFMLogManager();
    log_manager.writeOnLog("\nBegin Adam::findMin() procedure\n");

    // compute the current value
    double newf, newdf;
    this->_gradtargetfun->f(_x->getX(), newf, newdf);
    _x->setF(newf, newdf);
    _writeCurrentXInLog();

    //begin the minimization loop
    int step = 0;
    while ( _isNotConverged() )
        {
            log_manager.writeOnLog("\n\nAdam::findMin() Step " + std::to_string(step+1) + "\n");

            // compute the gradient
            this->_gradtargetfun->grad(_x->getX(), grad, graderr);
            _writeGradientInLog(grad, graderr);

            // compute the update
            for (int i=0; i<_ndim; ++i) {
                m[i] = _beta1 * m[i] + (1.-_beta1) * grad[i]; // Update biased first moment
                v[i] = _beta2 * v[i] + (1.-_beta2) * grad[i]*grad[i]; // Update biased second raw moment

                dx[i] = - afac * m[i] / (sqrt(v[i]) + _epsilon); // calculate updates
                _x->setX(i, _x->getX(i) + dx[i]); // update _x
            }

            _writeXUpdateInLog(dx);

            // compute the new value of the target function
            double newf, newdf;
            _gradtargetfun->f(_x->getX(), newf, newdf);
            _x->setF(newf, newdf);

            _writeCurrentXInLog();

            step++;
        }

    log_manager.writeOnLog("\nEnd Adam::findMin() procedure\n");
}
