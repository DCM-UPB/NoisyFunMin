#include "nfm/FIRE.hpp"

#include "nfm/LogManager.hpp"

#include <numeric>
#include <cmath>

namespace nfm
{

FIRE::FIRE(NoisyFunctionWithGradient * targetfun, const double dtinit, const double dtmax):
        NFM(targetfun), _dt0(std::max(0., std::min(dtmax, dtinit))), _dtmax(std::max(0., dtmax))
{
    if (!_flag_gradfun) {
        throw std::invalid_argument("[FIRE] FIRE optimization requires a target function with gradient.");
    }
    // overwrite defaults
    _flag_gradErrStop = false; // don't stop on noisy-low gradients, by default
}

// --- Minimization

void FIRE::_findMin()
{
    LogManager::logString("\nBegin FIRE::findMin() procedure\n");

    // helper variables
    const std::vector<double> &F = _grad.val; // convenience reference (the force)
    std::vector<double> v(_grad.size()); // velocity vector

    double dt = _dt0; // current time-step
    double alpha = _alpha0; // current mixing factor
    int Npos = 0; // number of steps since "grad.v" was negative

    //begin the minimization loop
    int iter = 0;
    while (true) {
        ++iter;
        if (LogManager::isLoggingOn()) { // else skip string construction
            LogManager::logString("\nFIRE::findMin() Step " + std::to_string(iter) + "\n");
        }

        // compute the gradient and current target
        this->_updateTarget();
        if (this->_shouldStop()) { break; } // we are done

        const double P = std::inner_product(F.begin(), F.end(), v.begin(), 0.); // P = F.v
        this->_mixVelocity(v, alpha);

        if (P > 0.) {
            if (++Npos > _Nmin) {
                dt = std::min(dt*_finc, _dtmax);
                alpha *= _falpha;
            }
        }
        else {
            Npos = 0;
            dt *= _fdec;
            alpha = _alpha0;
            std::fill(v.begin(), v.end(), 0.);
        }

        // make MD step to find the next position/velocity
        this->_doMDStep(v, dt);
    }

    /*if (_useAveraging) { // calculate the old value average as end result
        this->_averageOldValues(); // perform average and store it in last
    }*/
    LogManager::logString("\nEnd FIRE::findMin() procedure\n");
}

// --- Internal methods

void FIRE::_updateTarget()
{
    _last.f = _gradfun->fgrad(_last.x, _grad);
    this->_storeLastValue();
    this->_writeGradientToLog();
}

void FIRE::_mixVelocity(std::vector<double> &v, double alpha)
{
    const std::vector<double> &F = _grad.val; // we need only the values (the force)
    const double fnorm = sqrt(std::inner_product(F.begin(), F.end(), F.begin(), 0.));
    const double vnorm = sqrt(std::inner_product(v.begin(), v.end(), v.begin(), 0.));
    for (int i = 0; i < _ndim; ++i) {
        v[i] = (1. - alpha)*v[i] + alpha*vnorm*F[i]/fnorm;
    }
}

void FIRE::_doMDStep(std::vector<double> &v, double dt)
{
    // compute euler update
    for (int i = 0; i< _ndim; ++i) {
        v[i] += dt*_grad.val[i];
        _last.x[i] += dt*v[i];
    }
}
} // namespace nfm
