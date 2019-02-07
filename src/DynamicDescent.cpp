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

    stringstream s;
    s << "inertia: " << _inertia << endl;
    s << flush;
    log_manager.writeOnLog(s.str());
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
                _inertia = 0.;
                for (int i=0; i<this->_ndim; ++i)
                    {
                        _inertia+=grad[i]*grad[i];
                    }
                _inertia = _ndim / sqrt(_inertia);
            }

            // find the next position
            this->findNextX(grad);

            cont ++;
        }

    log_manager.writeNoisyValueInLog(_x, "Final position and target value");
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
    // compute the dot product between the normalized gradection and the old normalized gradection
    double old_new_direction_dot_product = 0;
    for (int i=0; i<_ndim; ++i){
        old_new_direction_dot_product += _old_norm_direction[i] * norm_grad[i];
    }
    // update the inertia and report it in the log
    _inertia = _inertia + 0.5 * _inertia * old_new_direction_dot_product;
    this->_writeInertiaInLog();

    // update _x
    double dx[_ndim];
    for (int i=0; i<_ndim; ++i){
        dx[i] = - _inertia * grad[i];
        _x->setX(i, _x->getX(i) + dx[i]);
    }
    this->_writeXUpdateInLog(dx);

    // store the grad for the next iteration
    for (int i=0; i<this->_ndim; ++i){
        _old_norm_direction[i] = norm_grad[i];
    }
}
