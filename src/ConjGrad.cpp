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

void ConjGrad::_findMin()
{
    using namespace std;

    LogManager::logString("\nBegin ConjGrad::findMin() procedure\n");

    //initialize the gradients
    std::vector<NoisyValue> gradh(static_cast<size_t>(_ndim));
    _last.f = _gradfun->fgrad(_last.x, gradh);
    this->_storeLastValue();
    this->_writeGradientToLog(gradh);

    if (this->_meaningfulGradient(gradh)) {
        // initialize gradient vector and old length
        std::vector<double> gradnew(static_cast<size_t>(_ndim)); // store negative gradient in gradnew
        double scalprodold = 0.;
        for (int i = 0; i < _ndim; ++i) {
            gradnew[i] = -gradh[i].value; // interested in following -gradient
            scalprodold += gradnew[i]*gradnew[i];
        }
        scalprodold = sqrt(scalprodold);

        // initialize the conjugate vectors
        std::vector<double> conjv(gradnew);
        this->_writeCGDirectionToLog(conjv, "Conjugated Vectors");

        //find new position
        this->_findNextX(conjv);

        //begin the minimization loop
        int cont = 0;
        while (!this->_shouldStop()) {
            ++cont;
            LogManager::logString("\n\nConjGrad::findMin() Step " + std::to_string(cont) + "\n");
            //cout << "x is in " << getX(0) << "   " << getX(1) << "   " << getX(2) << endl << endl;
            //evaluate the new gradient
            _gradfun->grad(_last.x, gradh);
            this->_writeGradientToLog(gradh);
            if (!this->_meaningfulGradient(gradh)) { break; }
            for (int i = 0; i < _ndim; ++i) { gradnew[i] = -gradh[i].value; } // negative gradient

            // compute the direction to follow for finding the next x
            //    if _use_conjgrad == true   ->   Conjugate Gradient
            //    else   ->   Steepest Descent
            switch (_cgmode) {
            case CGMode::RAW:
                // simply use as conjugate gradient the gradient (i.e. make a steepest descent!)
                for (int i = 0; i < _ndim; ++i) { std::copy(gradnew.begin(), gradnew.end(), conjv.begin()); }
                break;
            case CGMode::CGFR:
                //determine the new conjugate vector (with Fletcher-Reeves step)
                const double scalprodnew = std::inner_product(gradnew.begin(), gradnew.end(), gradnew.begin(), 0.);
                const double ratio = scalprodold != 0 ? scalprodnew/scalprodold : 0.;
                scalprodold = scalprodnew;
                for (int i = 0; i < _ndim; ++i) {
                    conjv[i] = gradnew[i] + conjv[i]*ratio;
                }
                break;
            }
            this->_writeCGDirectionToLog(conjv, "Conjugated vectors");

            //find new position
            this->_findNextX(conjv);
        }
    }

    LogManager::logString("\nEnd ConjGrad::findMin() procedure\n");
}


// --- Internal methods

void ConjGrad::_findNextX(const std::vector<double> &dir)
{
    // do line-minimization and store result in last
    _last = nfm::multiLineMin(*_targetfun, _last, dir, 0.25, 1.0, std::max(_epsx, m1d_default::XTOL), std::max(_epsf, m1d_default::FTOL)); // keep non-zero tol for line-search algos
    this->_storeLastValue();
}
} // namespace nfm