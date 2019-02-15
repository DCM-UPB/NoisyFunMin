#include "nfm/DynamicDescent.hpp"

#include "nfm/LogNFM.hpp"
#include "nfm/FunProjection1D.hpp"
#include "nfm/1DTools.hpp"

#include <iostream>
#include <sstream>
#include <cmath>
#include <stdexcept>


// --- Log

void DynamicDescent::_writeInertiaInLog(){
    using namespace std;

    NFMLogManager log_manager = NFMLogManager();
    log_manager.writeVectorInLog(_inertia, NULL, _ndim, 2, "Current inertia", "i");
}


// --- Minimization

void DynamicDescent::findMin(){
    using namespace std;

    NFMLogManager log_manager = NFMLogManager();
    log_manager.writeOnLog("\nBegin DynamicDescent::findMin() procedure\n");

    // clear old values
    this->_clearOldValues();

    //arrays to hold the gradients
    double grad[_ndim];
    double graderr[_ndim];

    //begin the minimization loop
    int cont = 0;
    while ( true )
        {

            // compute the gradient and current target
            double newf, newdf;
            this->_gradtargetfun->fgrad(_x->getX(), newf, newdf, grad, graderr);
            _x->setF(newf, newdf);

            if (this->_shouldStop(grad, graderr)) break;

            log_manager.writeOnLog("\n\nDynamicDescent::findMin() Step " + std::to_string(cont+1) + "\n");
            this->_writeCurrentXInLog();
            this->_writeGradientInLog(grad, graderr);

            // if it is the first iteration, initialise the inertia
            if (cont == 0){
                for (int i=0; i<_ndim; ++i) {
                    _inertia[i] = (grad[i]!=0) ? 1. / fabs(grad[i]) : 0.;
                }
            }

            // find the next position
            this->findNextX(grad);

            cont ++;
        }

    log_manager.writeNoisyValueInLog(_x, 1, "Final position and target value");
    log_manager.writeOnLog("\nEnd DynamicDescent::findMin() procedure\n");
}


// --- Internal methods

void DynamicDescent::findNextX(const double * grad)
{
    using namespace std;

    // compute the normalized gradection vector
    double norm_grad[_ndim];
    double sum = 0.;
    for (int i=0; i<_ndim; ++i){
        sum += grad[i]*grad[i];
    }
    sum = sqrt(sum);
    for (int i=0; i<_ndim; ++i){
        if (sum!=0.)
            norm_grad[i] = grad[i]/sum;
        else
            norm_grad[i] = 0.;
    }
    // update the inertia
    for (int i=0; i<_ndim; ++i){
        _inertia[i] += 0.5 * _inertia[i] * _old_norm_direction[i] * norm_grad[i];
    }
    // report it in the log
    this->_writeInertiaInLog();

    // update _x
    double dx[_ndim];
    for (int i=0; i<_ndim; ++i){
        dx[i] = - _step_size*_inertia[i] * grad[i];
        _x->setX(i, _x->getX(i) + dx[i]);
    }
    this->_writeXUpdateInLog(dx);

    // store the grad for the next iteration
    for (int i=0; i<this->_ndim; ++i){
        _old_norm_direction[i] = norm_grad[i];
    }
}
