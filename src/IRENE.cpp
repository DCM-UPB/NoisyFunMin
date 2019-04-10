#include "nfm/IRENE.hpp"

#include "nfm/LogManager.hpp"

#include <cmath>
#include <iostream>

namespace nfm
{

// --- Minimization

void IRENE::_findMin()
{
    LogManager::logString("\nBegin IRENE::findMin() procedure\n");

    // helper variables
    std::vector<double> v(_grad.size()); // velocity vector
    std::vector<double> ma(_grad.size()); // moving average acceleration
    NoisyGradient a(_ndim); // noisy mixed acceleration vector (used for MD)

    std::function<void()> update = [&]()
    { // MD force update lambda
        _last.f = _gradfun->fgrad(_last.x, _grad);
        md::computeAcceleration(_grad.val, _mi, a.val);
        md::computeAcceleration(_grad.err, _mi, a.err); // we need that later
        if (_beta > 0.) {
            for (int i = 0; i < _ndim; ++i) { // update averaged acceleration and mixed acceleration (for MD)
                ma[i] = _beta*ma[i] + (1. - _beta)*a.val[i];
                const double ISNR = std::min(1., a.err[i]/fabs(ma[i]));
                a.val[i] = (a.val[i] + ISNR*ISNR*ma[i])/(ISNR*ISNR + 1.); // mix raw and averaged gradient according to ISNR squared
                // note that we currently "assume" that the statistical error of a does not change due to this
            }
        }
    };
    md::MDView mdview{.x = _last.x, .v = v, .a = a.val, .update = update}; // references for MD integrator

    double dt = _dt0; // current time-step
    double alpha = _alpha0; // current mixing factor
    int Npos = 0; // number of steps since "F.v" was negative
    int Nmin = 0; // number of steps since dt = dtmin

    // initial step
    if (!this->_initializeMD(mdview, dt)) { return; } // return if shouldStop() already

    //begin the minimization loop
    int iter = 0;
    while (true) {
        ++iter;
        if (LogManager::isLoggingOn()) { // else skip string construction
            LogManager::logString("\nIRENE::findMin() Step " + std::to_string(iter) + "\n");
        }

        // compute the gradient and current target
        this->_updateTarget(mdview, dt);
        std::transform(_grad.err.begin(), _grad.err.end(), _mi.begin(), a.err.begin(), std::multiplies<>()); // update a.err
        if (this->_isNDtMinReached(Nmin) || this->_shouldStop()) { break; } // we are done

        // compute P, which is a NoisyValue in this algorithm
        NoisyValue P{};
        P.val = std::inner_product(v.begin(), v.end(), a.val.begin(), 0.); // P = v.a
        // error estimation based on error propagation
        for (int i = 0; i < _ndim; ++i) {
            P.err += pow(v[i]*a.err[i], 2); // dP/da * a_err
        }
        P.err = sqrt(P.err);

        // acceleration / velocity mixing
        const double vnorm = sqrt(std::inner_product(v.begin(), v.end(), v.begin(), 0.));
        const double anorm = sqrt(std::inner_product(a.val.begin(), a.val.end(), a.val.begin(), 0.));
        for (int i = 0; i < _ndim; ++i) {
            v[i] = (1. - alpha)*v[i] + alpha*vnorm*a.val[i]/anorm;
        }

        // check P (using noisy comparison)
        if (P > 0.) { // we are most likely going downhill
            if (++Npos > _Nwait) { // then increase dt
                dt = std::min(dt*_finc, _dtmax);
                Nmin = 0; // we have increased dt
                alpha *= _falpha;
            }
        }
        else if (P < 0.) { // we are most likely going uphill
            Npos = 0;
            dt = std::max(dt*_fdec, _dtmin);
            alpha = _alpha0;

            if (_flag_fullFreeze) { // freeze the system completely
                std::fill(v.begin(), v.end(), 0.);
            }
            else { // freeze selectively
                for (int i = 0; i < _ndim; ++i) {
                    if (a[i]*v[i] < 0.) { // using two noisy overloads here!
                        v[i] = 0.;
                    }
                }
            }
        }

        if (dt == _dtmin) { ++Nmin; }
    }

    LogManager::logString("\nEnd IRENE::findMin() procedure\n");
}
} // namespace nfm