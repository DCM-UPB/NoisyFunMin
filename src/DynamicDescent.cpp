#include "DynamicDescent.hpp"

#include "LogNFM.hpp"
#include "FunProjection1D.hpp"
#include "1DTools.hpp"

#include <iostream>
#include <sstream>
#include <cmath>
#include <stdexcept>


// --- Log

void DynamicDescent::_writeInertiaInLog(){
    using namespace std;

    NFMLogManager log_manager = NFMLogManager();

    stringstream s;
    s << endl << "inertia: " << _inertia << endl;
    s << flush;
    log_manager.writeOnLog(s.str());
}


// --- Minimization

void DynamicDescent::findMin(){
    using namespace std;

    NFMLogManager log_manager = NFMLogManager();
    log_manager.writeOnLog("\nBegin DynamicDescent::findMin() procedure\n");

    //initialize the gradients
    double grad[_ndim];
    double graderr[_ndim];

    // compute the current value
    double newf, newdf;
    this->_gradtargetfun->f(_x->getX(), newf, newdf);
    _x->setF(newf, newdf);
    this->_writeCurrentXInLog();


    //begin the minimization loop
    int cont = 0;
    while ( this->_isNotConverged() )
        {
            // compute the gradient
            this->_gradtargetfun->grad(_x->getX(), grad, graderr);
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
            this->_writeCurrentXInLog();

            cont ++;
        }

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

    // compute the new value of the target function
    double newf, newdf;
    _gradtargetfun->f(_x->getX(), newf, newdf);

    // store the new function value in _x
    _x->setF(newf, newdf);

    // store the grad for the next iteration
    for (int i=0; i<this->_ndim; ++i){
        _old_norm_direction[i] = norm_grad[i];
    }
}
