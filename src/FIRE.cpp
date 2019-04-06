#include "nfm/FIRE.hpp"

#include "nfm/LogManager.hpp"

#include <numeric>
#include <cmath>

#include <iostream>
namespace nfm
{

FIRE::FIRE(NoisyFunctionWithGradient * targetfun, const double dtmax, const double dt0):
        NFM(targetfun), _dtmax(std::max(0., dtmax)), _dt0((dt0 > 0.) ? std::min(dt0, dtmax) : 0.1*dtmax)
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
    std::vector<double> v(F.size()); // velocity vector

    double dt = _dt0; // current time-step
    double alpha = _alpha0; // current mixing factor
    int Npos = 0; // number of steps since "F.v" was negative

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

        // update vector lengths and P
        const double vnorm = sqrt(std::inner_product(v.begin(), v.end(), v.begin(), 0.));
        const double fnorm = sqrt(std::inner_product(F.begin(), F.end(), F.begin(), 0.));
        const double P = std::inner_product(F.begin(), F.end(), v.begin(), 0.);

        // modify velocity
        for (int i = 0; i < _ndim; ++i) {
            v[i] = (1. - alpha)*v[i] + alpha*vnorm*F[i]/fnorm;
        }

        // check P
        if (P > 0.) { // we are going downhill
            if (++Npos > _Nmin) { // then increase dt
                dt = std::min(dt*_finc, _dtmax);
                alpha *= _falpha;
            }
        }
        else { // we are going uphill
            Npos = 0;
            dt *= _fdec;
            alpha = _alpha0;
            std::fill(v.begin(), v.end(), 0.); // freeze the system
        }

        // make MD step to find the next position/velocity
        this->_doMDStep(v, dt);
    }

    LogManager::logString("\nEnd FIRE::findMin() procedure\n");
}

// --- Internal methods

void FIRE::_updateTarget()
{
    _last.f = _gradfun->fgrad(_last.x, _grad);
    this->_storeLastValue();
    this->_writeGradientToLog();
}

void FIRE::_doMDStep(std::vector<double> &v, double dt)
{
    // compute explicit Euler update
    for (int i = 0; i < _ndim; ++i) {
        v[i] += dt*_grad.val[i];
        _last.x[i] += dt*v[i];
    }
}
} // namespace nfm
