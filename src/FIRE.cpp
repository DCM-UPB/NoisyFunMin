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
    std::vector<double> v(F.size()); // velocity vector
    std::vector<double> m(F.size()); // gradient running average (only for beta>0)

    double dt = _dt0; // current time-step
    double alpha = _alpha0; // current mixing factor
    int Npos = 0; // number of steps since "m.v" was negative

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

        // update vector lengths
        const double vnorm = sqrt(std::inner_product(v.begin(), v.end(), v.begin(), 0.));
        const double fnorm = sqrt(std::inner_product(F.begin(), F.end(), F.begin(), 0.));

        // calculate scalar product P = F.v or P = m.v (if beta > 0)
        double P;
        if (_beta > 0) {
            for (int i = 0; i < _ndim; ++i) {
                m[i] = _beta*m[i] + (1. - _beta)*F[i];
            }
            const double mnorm = sqrt(std::inner_product(m.begin(), m.end(), m.begin(), 0.));
            P = std::inner_product(m.begin(), m.end(), v.begin(), 0.);
            P = (P != 0.) ? P/(mnorm*vnorm) : 0.; // P would be 0 if mnorm or vnorm were 0
        }
        else {
            P = std::inner_product(F.begin(), F.end(), v.begin(), 0.);
            P = (P != 0.) ? P/(fnorm*vnorm) : 0.; // P would be 0 if fnorm or vnorm were 0
        }

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
            if (_flag_soft) { // freeze softly
                for (auto &vi : v) {
                    vi *= (1. - fabs(P));
                }
            }
            else { // freeze completely (default)
                std::fill(v.begin(), v.end(), 0.);
            }
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
