#include "nfm/ConjGrad.hpp"

#include "nfm/1DTools.hpp"
#include "nfm/LogManager.hpp"

#include <cmath>
#include <iostream>
#include <numeric>

namespace nfm
{

// --- Logging

void ConjGrad::_writeCGDirectionToLog(const std::vector<double> &dir, const std::string &name) const
{
    LogManager::logVector(dir, LogLevel::VERBOSE, name, "g");
}


// --- Minimization

void ConjGrad::findMin()
{
    using namespace std;

    LogManager::logString("\nBegin ConjGrad::findMin() procedure\n");

    // clear old values
    this->_clearOldValues();

    //initialize the gradients
    std::vector<NoisyValue> gradh(static_cast<size_t>(_ndim));
    _last.f = _gradfun->fgrad(_last.x, gradh);
    this->_storeOldValue();

    // store inverted gradient in gradold
    std::vector<double> gradold(static_cast<size_t>(_ndim));
    for (int i = 0; i<_ndim; ++i) { gradold[i] = -gradh[i].value; } //interested in following -gradient

    this->_writeCurrentXToLog();
    this->_writeGradientToLog(gradh);

    if (this->_meaningfulGradient(gradh)) {
        // internal vectors
        std::vector<double> gradnew(static_cast<size_t>(_ndim)); // will be used to store new grad
        std::vector<double> conjv(gradold); //initialize the conjugate vectors
        this->_writeCGDirectionToLog(conjv, "Conjugated Vectors");

        //find new position
        double deltatargetfun, deltax;
        this->findNextX(conjv, deltatargetfun, deltax);
        this->_writeCurrentXToLog();

        //begin the minimization loop
        //cout << "deltatargetfunction = " << deltatargetfun << "   " << _epsf << endl;
        //cout << "deltax = " << deltax << "   " << _epsx << endl << endl;
        int cont = 0;
        while ((deltatargetfun >= _epsf) && (deltax >= _epsx)) {
            LogManager::logString("\n\nConjGrad::findMin() Step " + std::to_string(cont + 1) + "\n");
            //cout << "x is in " << getX(0) << "   " << getX(1) << "   " << getX(2) << endl << endl;
            //evaluate the new gradient (and function value)
            _last.f = _gradfun->fgrad(_last.x, gradh);
            this->_writeGradientToLog(gradh);
            if (!this->_meaningfulGradient(gradh)) { break; }
            for (int i = 0; i<_ndim; ++i) { gradnew[i] = -gradh[i].value; } // negative gradient

            // compute the direction to follow for finding the next x
            //    if _use_conjgrad == true   ->   Conjugate Gradient
            //    else   ->   Steepest Descent
            if (_use_conjgrad) {
                //determine the new conjugate vector
                const double scalprodnew = std::inner_product(gradnew.begin(), gradnew.end(), gradnew.begin(), 0.);
                const double scalprodold = std::inner_product(gradold.begin(), gradold.end(), gradold.begin(), 0.);
                const double ratio = scalprodnew/scalprodold;
                for (int i = 0; i < _ndim; ++i) { conjv[i] = gradnew[i] + conjv[i]*ratio; }
            }
            else {
                // simply use as conjugate gradient the gradient (i.e. make a steepest descent!)
                for (int i = 0; i< _ndim; ++i) { std::copy(gradnew.begin(), gradnew.end(), conjv.begin()); }
            }
            this->_writeCGDirectionToLog(conjv, "Conjugated vectors");

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

    LogManager::logString("\nEnd ConjGrad::findMin() procedure\n");
}


// --- Internal methods

void ConjGrad::findNextX(const std::vector<double> &dir, double &deltatargetfun, double &deltax)
{
    using namespace std;

    NoisyIOPair old = _last; // on the last gradient calculation f was stored as well
    _last = nfm::multiLineMinimization(*_targetfun, _last.x, dir, _epsf); // store line-minimization result in last
    //compute the two deltas
    deltatargetfun = std::max(0., fabs(_last.f.value - old.f.value) - _last.f.error - old.f.error);
    std::transform (old.x.begin(), old.x.end(), _last.x.begin(), old.x.begin(), std::minus<>()); // old.x = old.x-last.x
    deltax = std::inner_product(old.x.begin(), old.x.end(), old.x.begin(), 0.); // delta = old.x . old.x
    deltax = sqrt(deltax); // distance from (original) old.x to last.x
    //cout << "x is in " << this->getX(0) << "   " << this->getX(1) << "   " << this->getX(2) << endl;
}
} // namespace nfm