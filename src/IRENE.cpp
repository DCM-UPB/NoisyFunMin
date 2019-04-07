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
    //std::vector<double> ma(_grad.size()); // moving average acceleration
    NoisyGradient a(_ndim); // noisy acceleration vector (F*mi)
    md::MDView mdview{.mi = _mi, .x = _last.x, .v = v, .a = a.val, .F = _grad}; // references for MD integrator

    double dt = _dt0; // current time-step
    double alpha = _alpha0; // current mixing factor
    int Npos = 0; // number of steps since "F.v" was negative
    int Nmin = 0; // number of steps since dt = dtmin

    // initial step (velocities are kept 0)
    if (!this->_initializeMD(mdview, dt)) { return; } // return if shouldStop() already

    //begin the minimization loop
    int iter = 0;
    while (true) { // a better while (True)
        ++iter;
        if (LogManager::isLoggingOn()) { // else skip string construction
            LogManager::logString("\nIRENE::findMin() Step " + std::to_string(iter) + "\n");
        }
        //std::cout << std::endl << "IRENE step " << iter << ":" << std::endl;

        // compute the gradient and current target
        this->_updateTarget(mdview, dt);
        std::transform(_grad.err.begin(), _grad.err.end(), _mi.begin(), a.err.begin(), std::multiplies<>()); // update a.err
        if (this->_isNDtMinReached(Nmin) || this->_shouldStop()) { break; } // we are done

        // update vector lengths and P
        const double vnorm = sqrt(std::inner_product(v.begin(), v.end(), v.begin(), 0.));
        const double anorm = sqrt(std::inner_product(a.val.begin(), a.val.end(), a.val.begin(), 0.));

        // compute P, which is a NoisyValue in this algorithm
        NoisyValue P{};
        P.val = std::inner_product(v.begin(), v.end(), a.val.begin(), 0.); // P = v.a
        // error estimation based on error propagation
        for (int i = 0; i < _ndim; ++i) {
            P.err += pow(v[i]*a.err[i], 2); // dP/da * a_err
        }
        P.err = sqrt(P.err);

        //std::cout << "vnorm " << vnorm << ", anorm " << anorm << ", P" << P << std::endl;

        // velocity mixing
        for (int i = 0; i < _ndim; ++i) {
            v[i] = (1. - alpha)*v[i] + alpha*vnorm*a.val[i]/anorm;
        }

        // check P (using noisy comparison)
        if (P > 0.) { // we are definitely going downhill
            if (++Npos > _Nwait) { // then increase dt
                dt = std::min(dt*_finc, _dtmax);
                Nmin = 0; // we have increased dt
                alpha *= _falpha;
            }
        }
        else if (P < 0.) { // we are definitely going uphill
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
        //std::cout << "dt " << dt << std::endl;
    }

    LogManager::logString("\nEnd IRENE::findMin() procedure\n");
}
    /*
        /// we use the averaged gradient
            for (int i = 0; i < _ndim; ++i) {
                m.set(i, m[i]*_beta + _grad[i]*(1. - _beta));
            }
            P.val = std::inner_product(m.val.begin(), m.val.end(), v.begin(), 0.);
            for (int i = 0; i < _ndim; ++i) {
                P.err += pow(m.err[i]*v[i], 2);
            }

                for (int i=0; i<_ndim; ++i) {
                    const double p_err = fabs(NoisyValue::getSigmaLevel()*Ferr[i]*v[i]); // err on Di*vi
                    if (F[i]*v[i] < -p_err) {
                        const double vold = v[i];
                        //v[i] = 0.;
                        v[i] = (F[i] != 0.) ? v[i] = p_err/F[i] : 0.;
                        //std::cout << "flip v" << i << ": " << vold << " -> " << v[i] << std::endl;
                    }
   */
} // namespace nfm