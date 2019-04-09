#include "nfm/ConjGrad.hpp"

#include "nfm/LogManager.hpp"

#include <cmath>
#include <numeric>

namespace nfm
{

// --- Constructor

ConjGrad::ConjGrad(NoisyFunctionWithGradient * targetfun, const CGMode cgmode, const MLMParams params):
        NFM(targetfun), _cgmode(cgmode), _mlmParams(params)
{
    if (!_flag_gradfun) {
        throw std::invalid_argument("[ConjGrad] Conjugate Gradient optimization requires a target function with gradient.");
    }
    // overwrite defaults
    this->setMaxNConstValues(1); // don't use the check by default
    this->setEpsX(m1d_detail::STD_XTOL); // because this means we stop on reject line search
    this->setEpsF(m1d_detail::STD_FTOL); // and this means we stop if it didn't improve target significantly (beyond tol+errors)
}

// --- Logging

void ConjGrad::_writeCGDirectionToLog(const std::vector<double> &dir, const std::string &name) const
{
    LogManager::logVector(dir, LogLevel::VERBOSE, name, "g");
}

// --- Minimization

void ConjGrad::_findMin()
{
    // --- Starting Position

    LogManager::logString("\nBegin ConjGrad::findMin() procedure\n");

    // obtain the initial function value and gradient (uphill)
    bool flag_cont = this->_computeGradient(true); // compute grad and value
    if (!flag_cont) { return; } // return early


    // --- Initialize CG

    // initialize gradient vectors and length
    std::vector<double> &gradnew = _grad.val; // store reference to gradient values
    std::vector<double> conjv = gradnew; // stores the conjugate vectors, initialize with raw gradient
    std::vector<double> gradold; // the previous inverted gradients (only used for Polak-Ribiere CG)

    // save old gradient for PR-CG
    if (_cgmode == CGMode::CGPR || _cgmode == CGMode::CGPR0) {
        gradold = gradnew; // initialize old gradient
    }
    // the denominator of CG update ratio
    double gdot_old = std::inner_product(gradnew.begin(), gradnew.end(), gradnew.begin(), 0.);

    // find initial new position
    LogManager::logString("\nConjGrad::findMin() Step 1\n");
    this->_findNextX(conjv);


    // --- Main CG Loop

    int iter = 1;
    while (!this->_shouldStop()) {
        ++iter;
        if (LogManager::isLoggingOn()) { // else skip string construction
            LogManager::logString("\nConjGrad::findMin() Step " + std::to_string(iter) + "\n");
        }

        // evaluate the new gradient
        flag_cont = this->_computeGradient(false);
        if (!flag_cont) { return; } // gradient is only noise

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
                for (int i = 0; i < _ndim; ++i) { prprod += gradnew[i]*(gradnew[i] - gradold[i]); }
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

bool ConjGrad::_computeGradient(const bool flag_value)
{
    if (flag_value) { // value and gradient
        _last.f = _gradfun->fgrad(_last.x, _grad);
        this->_storeLastValue();
    }
    else { // only gradient
        _gradfun->grad(_last.x, _grad);
    }
    this->_writeGradientToLog();
    if (this->_isGradNoisySmall()) { // we directly check and print the exit message here
        LogManager::logString("\nEnd ConjGrad::findMin() procedure\n");
        return false;
    }
    return true;
}

void ConjGrad::_findNextX(const std::vector<double> &dir)
{
    // use NFM tolerances for MLM
    _mlmParams.epsx = this->getEpsX();
    _mlmParams.epsf = this->getEpsF();

    // do line-minimization and store result in last
    _last = nfm::multiLineMin(*_targetfun, _last, dir, _mlmParams);
    this->_storeLastValue();
}
} // namespace nfm
