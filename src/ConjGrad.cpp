#include "nfm/ConjGrad.hpp"

#include "nfm/1DTools.hpp"
#include "nfm/LogNFM.hpp"

#include <cmath>
#include <iostream>
#include <numeric>

namespace nfm
{

// --- Logging

void ConjGrad::_writeCGDirectionInLog(const double * dir, const std::string &name)
{
    NFMLogManager log_manager;
    log_manager.writeVectorInLog(dir, nullptr, _ndim, 2, name, "g");
}


// --- Minimization

void ConjGrad::findMin()
{
    using namespace std;

    NFMLogManager log_manager = NFMLogManager();
    log_manager.writeOnLog("\nBegin ConjGrad::findMin() procedure\n");

    // clear old values
    this->_clearOldValues();

    //initialize the gradients
    double gradold[_ndim];
    double graderr[_ndim];
    double newf, newerr;

    this->_gradtargetfun->fgrad(_x->getX(), newf, newerr, gradold, graderr);
    _x->setF(newf, newerr);
    this->_storeOldValue();

    this->_writeCurrentXInLog();
    this->_writeGradientInLog(gradold, graderr);

    if (this->_meaningfulGradient(gradold, graderr)) {
        double gradnew[_ndim];
        double conjv[_ndim]; //interested in following -gradient
        for (int i = 0; i < _ndim; ++i) { gradold[i] = -gradold[i]; }
        //inizialize the conjugate vectors
        for (int i = 0; i < _ndim; ++i) { conjv[i] = gradold[i]; }
        this->_writeCGDirectionInLog(conjv, "Conjugated Vectors");
        //find new position
        double deltatargetfun, deltax;
        this->findNextX(conjv, deltatargetfun, deltax);

        this->_writeCurrentXInLog();

        //begin the minimization loop
        //cout << "deltatargetfunction = " << deltatargetfun << "   " << _epstargetfun << endl;
        //cout << "deltax = " << deltax << "   " << _epsx << endl << endl;
        int cont = 0;
        while ((deltatargetfun >= _epstargetfun) && (deltax >= _epsx)) {
            log_manager.writeOnLog("\n\nConjGrad::findMin() Step " + std::to_string(cont + 1) + "\n");
            //cout << "x is in " << getX(0) << "   " << getX(1) << "   " << getX(2) << endl << endl;
            //evaluate the new gradient
            this->_gradtargetfun->grad(_x->getX(), gradnew, graderr);
            this->_writeGradientInLog(gradnew, graderr);
            for (int i = 0; i < _ndim; ++i) { gradnew[i] = -gradnew[i]; }
            if (!this->_meaningfulGradient(gradnew, graderr)) { break; }
            // compute the direction to follow for finding the next x
            //    if _use_conjgrad == true   ->   Conjugate Gradient
            //    else   ->   Steepest Descent
            if (_use_conjgrad) {
                //determine the new conjugate vector
                const double scalprodnew = std::inner_product(gradnew, gradnew + _ndim, gradnew, 0.);
                const double scalprodold = std::inner_product(gradold, gradold + _ndim, gradold, 0.);
                const double ratio = scalprodnew/scalprodold;
                for (int i = 0; i < _ndim; ++i) { conjv[i] = gradnew[i] + conjv[i]*ratio; }
            }
            else {
                // simply use as conjugate gradient the gradient (i.e. make a steepest descent!)
                std::copy(gradnew, gradnew + _ndim, conjv);
            }
            this->_writeCGDirectionInLog(conjv, "Conjugated vectors");
            //find new position
            this->findNextX(conjv, deltatargetfun, deltax);
            //cout << "deltatargetfunction = " << deltatargetfun << "   " << _epstargetfun << endl;
            //cout << "deltax = " << deltax << "   " << _epsx << endl << endl;

            this->_writeCurrentXInLog();

            _storeOldValue();
            if (this->_isConverged()) { break; }
            cont++;
        }
    }

    log_manager.writeOnLog("\nEnd ConjGrad::findMin() procedure\n");
}


// --- Internal methods

void ConjGrad::findNextX(const double * dir, double &deltatargetfun, double &deltax)
{
    using namespace std;

    //project the original multidimensional wave function into a one-dimensional function
    auto * proj1d = new FunProjection1D(_ndim, _x->getX(), dir, _targetfun);
    //determine the initial bracket
    NoisyValue a(1), b(1), c(1);
    a.setX(0.);
    double newf, dnewf;
    proj1d->f(a.getX(), newf, dnewf);
    a.setF(newf, dnewf);
    nfm::findBracket(proj1d, a, b, c);
    //find the minimum in the bracket
    nfm::parabgoldMinimization(proj1d, _epstargetfun, a, b, c);
    //get the x corresponding to the found b
    proj1d->getVecFromX(b.getX(0), _x->getX());
    _x->setF(b.getF(), b.getDf());
    //compute the two deltas
    deltatargetfun = fabs(b.getF() - newf) - dnewf - b.getDf();
    deltax = std::inner_product(dir, dir + _ndim, dir, 0.);
    deltax = fabs(b.getX(0)*sqrt(deltax)); //*NORMA VETTORE dir
    //cout << "x is in " << this->getX(0) << "   " << this->getX(1) << "   " << this->getX(2) << endl;
    delete proj1d;
}
} // namespace nfm