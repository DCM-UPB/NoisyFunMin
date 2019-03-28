#include "nfm/ConjGrad.hpp"

#include "nfm/1DTools.hpp"
#include "nfm/LogManager.hpp"

#include <cmath>
#include <iostream>
#include <numeric>

namespace nfm
{

// --- Constructor

ConjGrad::ConjGrad(NoisyFunctionWithGradient * targetfun, const int max_n_const_values,
                   const double stepSize, const int max_n_bracketing):
        NFM(targetfun, max_n_const_values), _stepSize(stepSize),
        _max_n_bracketing(std::max(1, max_n_bracketing)), _cgmode(CGMode::CGFR) {}

// --- Logging

void ConjGrad::_writeCGDirectionToLog(const std::vector<double> &dir, const std::string &name) const
{
    LogManager::logVector(dir, LogLevel::VERBOSE, name, "g");
}

// --- Minimization

void ConjGrad::_findMin()
{
    using namespace std;

    // --- Starting Position

    LogManager::logString("\nBegin ConjGrad::findMin() procedure\n");

    // obtain the initial function value and gradient (uphill)
    std::vector<NoisyValue> gradh(static_cast<size_t>(_ndim));
    _last.f = _gradfun->fgrad(_last.x, gradh);
    this->_storeLastValue();
    this->_writeGradientToLog(gradh);

    // initial sanity check
    if (!this->_meaningfulGradient(gradh)) {
        LogManager::logString("\nEnd ConjGrad::findMin() procedure\n");
        return;
    }

    // --- Initialize CG

    // initialize gradient vectors and length
    std::vector<double> gradnew(gradh.size()); // stores new raw gradients with inverted sign
    std::vector<double> conjv(gradh.size()); // stores the conjugate vectors
    std::vector<double> gradold; // the previous inverted gradients (only used for Polak-Ribiere CG)

    for (int i = 0; i < _ndim; ++i) {
        gradnew[i] = -gradh[i].value; // interested in following -gradient
        conjv[i] = gradnew[i]; // initialize CG vector with initial gradient
    }
    // save old gradient for PR-CG
    if (_cgmode == CGMode::CGPR || _cgmode == CGMode::CGPR0) {
        gradold = gradnew; // initialize old gradient
    }
    // the denominator of CG update ratio
    double gdot_old = std::inner_product(gradnew.begin(), gradnew.end(), gradnew.begin(), 0.);

    // find initial new position
    this->_findNextX(conjv);


    // --- Main CG Loop

    int cont = 0;
    while (!this->_shouldStop()) {
        ++cont;
        LogManager::logString("\n\nConjGrad::findMin() Step " + std::to_string(cont) + "\n");

        // evaluate the new gradient
        _gradfun->grad(_last.x, gradh);
        this->_writeGradientToLog(gradh);
        if (!this->_meaningfulGradient(gradh)) { break; }
        for (int i = 0; i < _ndim; ++i) { gradnew[i] = -gradh[i].value; } // negative gradient

        // compute the next direction to follow
        if (_cgmode == CGMode::NOCG) { // use raw gradient (i.e. steepest descent)
            std::copy(gradnew.begin(), gradnew.end(), conjv.begin());
        }
        else { // use conjugate gradients
            const double gdot_new = std::inner_product(gradnew.begin(), gradnew.end(), gradnew.begin(), 0.);

            double ratio; // CG update factor
            if (_cgmode == CGMode::CGFR) { // Fletcher-Reeves CG
                ratio = gdot_old != 0 ? gdot_new/gdot_old : 0.;
            }
            else { // Polak-Ribiere CG
                double prprod = 0.;
                for (int i = 0; i < _ndim; ++i) { prprod += pow(gradnew[i] - gradold[i], 2); }
                ratio = gdot_old != 0 ? prprod/gdot_old : 0.;
                if (_cgmode == CGMode::CGPR0) { ratio = std::max(0., ratio); } // CG reset
                gradold = gradnew; // gradient old to new
            }
            gdot_old = gdot_new; // gdot old to new

            // update conjugate gradients
            for (int i = 0; i < _ndim; ++i) {
                conjv[i] = gradnew[i] + ratio*conjv[i];
            }
            this->_writeCGDirectionToLog(conjv, "Conjugated vectors");
        }

        // find new position and continue loop
        this->_findNextX(conjv);
    }

    LogManager::logString("\nEnd ConjGrad::findMin() procedure\n");
}


// --- Internal methods

void ConjGrad::_findNextX(const std::vector<double> &dir)
{
    // do line-minimization and store result in last
    _last = nfm::multiLineMin(*_targetfun, _last, dir, 0.5*_stepSize /*backstep*/, _stepSize/*fwdstep*/, _max_n_bracketing,
                              std::max(_epsx, m1d_default::XTOL), std::max(_epsf, m1d_default::FTOL)); // keep non-zero tol for line-search algos
    this->_storeLastValue();
}
} // namespace nfm