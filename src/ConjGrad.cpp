#include "nfm/ConjGrad.hpp"

#include "nfm/1DTools.hpp"
#include "nfm/LogManager.hpp"

#include <cmath>
#include <iostream>
#include <numeric>

namespace nfm
{

// --- Logging

void ConjGrad::_writeCGDirectionInLog(const double * dir, const std::string &name)
{
    LogManager log_manager;
    log_manager.writeVectorInLog(dir, nullptr, _ndim, 2, name, "g");
}


// --- Minimization

void ConjGrad::findMin()
{
    using namespace std;

    LogManager log_manager = LogManager();
    log_manager.logString("\nBegin ConjGrad::findMin() procedure\n");

    // clear old values
    this->_clearOldValues();

    //initialize the gradients
    double gradold[_ndim];
    double graderr[_ndim];
    double newf, newerr;

    this->_gradfun->fgrad(_last->getX(), newf, newerr, gradold, graderr);
    _last->setF(newf, newerr);
    this->_storeOldValue();

    this->_writeCurrentXToLog();
    this->_writeGradientToLog(gradold, graderr);

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

        this->_writeCurrentXToLog();

        //begin the minimization loop
        //cout << "deltatargetfunction = " << deltatargetfun << "   " << _epsf << endl;
        //cout << "deltax = " << deltax << "   " << _epsx << endl << endl;
        int cont = 0;
        while ((deltatargetfun >= _epsf) && (deltax >= _epsx)) {
            log_manager.logString("\n\nConjGrad::findMin() Step " + std::to_string(cont + 1) + "\n");
            //cout << "x is in " << getX(0) << "   " << getX(1) << "   " << getX(2) << endl << endl;
            //evaluate the new gradient
            this->_gradfun->grad(_last->getX(), gradnew, graderr);
            this->_writeGradientToLog(gradnew, graderr);
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
            //cout << "deltatargetfunction = " << deltatargetfun << "   " << _epsf << endl;
            //cout << "deltax = " << deltax << "   " << _epsx << endl << endl;

            this->_writeCurrentXToLog();

            _storeOldValue();
            if (this->_isConverged()) { break; }
            cont++;
        }
    }

    log_manager.logString("\nEnd ConjGrad::findMin() procedure\n");
}


// --- Internal methods

void ConjGrad::findNextX(const double * dir, double &deltatargetfun, double &deltax)
{
    using namespace std;

    //project the original multidimensional wave function into a one-dimensional function
    auto * proj1d = new FunProjection1D(_ndim, _last->getX(), dir, _targetfun);
    //determine the initial bracket
    NoisyValue a(1), b(1), c(1);
    a.setX(0.);
    double newf, dnewf;
    proj1d->f(a.getX(), newf, dnewf);
    a.setF(newf, dnewf);
    nfm::findBracket(proj1d, a, b, c);
    //find the minimum in the bracket
    nfm::parabgoldMinimization(proj1d, _epsf, a, b, c);
    //get the x corresponding to the found b
    proj1d->getVecFromX(b.getX(0), _last->getX());
    _last->setF(b.getF(), b.getDf());
    //compute the two deltas
    deltatargetfun = fabs(b.getF() - newf) - dnewf - b.getDf();
    deltax = std::inner_product(dir, dir + _ndim, dir, 0.);
    deltax = fabs(b.getX(0)*sqrt(deltax)); //*NORMA VETTORE dir
    //cout << "x is in " << this->getX(0) << "   " << this->getX(1) << "   " << this->getX(2) << endl;
    delete proj1d;
}
} // namespace nfm