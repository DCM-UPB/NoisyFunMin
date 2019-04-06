#include "nfm/FIRE.hpp"

#include "nfm/LogManager.hpp"

#include <cmath>


namespace nfm
{

FIRE::FIRE(NoisyFunctionWithGradient * targetfun, const double dtmax, const double dt0):
        NFM(targetfun), _dtmax(std::max(0., dtmax)), _dt0((dt0 > 0.) ? std::min(dt0, dtmax) : 0.1*dtmax)
{
    if (!_flag_gradfun) {
        throw std::invalid_argument("[FIRE] FIRE optimization requires a target function with gradient.");
    }
    _mi.assign(_grad.size(), 1.); // inverse masses default to 1
    // overwrite defaults
    _flag_gradErrStop = false; // don't stop on noisy-low gradients, by default
}

// --- Minimization

void FIRE::_findMin()
{
    LogManager::logString("\nBegin FIRE::findMin() procedure\n");

    // helper variables
    std::vector<double> v(_grad.size()); // velocity vector
    std::vector<double> a(_grad.size()); // acceleration vector (F*mi)
    md::MDView mdview{.mi = _mi, .x = _last.x, .v = v, .a = a, .F = _grad}; // references for MD integrator

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
            LogManager::logString("\nFIRE::findMin() Step " + std::to_string(iter) + "\n");
        }

        // compute the gradient and current target
        this->_updateTarget(mdview, dt);
        if ((_Ndtmin > 0 && Nmin > _Ndtmin) || this->_shouldStop()) { break; } // we are done

        // update vector lengths and P
        const double vnorm = sqrt(std::inner_product(v.begin(), v.end(), v.begin(), 0.));
        const double anorm = sqrt(std::inner_product(a.begin(), a.end(), a.begin(), 0.));
        const double P = std::inner_product(a.begin(), a.end(), v.begin(), 0.);

        // velocity mixing
        for (int i = 0; i < _ndim; ++i) {
            v[i] = (1. - alpha)*v[i] + alpha*vnorm*a[i]/anorm;
        }

        // check P
        if (P > 0.) { // we are going downhill
            if (++Npos > _Nwait) { // then increase dt
                dt = std::min(dt*_finc, _dtmax);
                Nmin = 0; // we have increased dt
                alpha *= _falpha;
            }
        }
        else { // we are going uphill
            Npos = 0;
            dt = std::max(dt*_fdec, _dtmin);
            if (dt == _dtmin) { ++Nmin; }
            alpha = _alpha0;

            if (_flag_fullFreeze) { // freeze the system completely
                std::fill(v.begin(), v.end(), 0.);
            }
            else { // freeze selectively
                for (int i = 0; i < _ndim; ++i) {
                    if (a[i]*v[i] < 0.) {
                        v[i] = 0.;
                    }
                }
            }
        }
    }

    LogManager::logString("\nEnd FIRE::findMin() procedure\n");
}

// --- Internal methods

bool FIRE::_initializeMD(md::MDView &view, const double dt)
{
    LogManager::logString("\nFIRE::findMin() Initial Step\n");
    _last.f = _gradfun->fgrad(_last.x, _grad);
    md::computeAcceleration(view); // store initial acceleration
    this->_storeLastValue();
    this->_writeGradientToLog();
    if (this->_shouldStop()) { // we print termination message already
        LogManager::logString("\nEnd FIRE::findMin() procedure\n");
        return false;
    }
    return true; // we can start the algorithm
}

void FIRE::_updateTarget(md::MDView &view, const double dt)
{
    _last.f = md::doMDStep(_mdi, *_gradfun, view, dt);
    this->_storeLastValue();
    this->_writeGradientToLog();
}

// --- Public

void FIRE::setMasses(const std::vector<double> &m)
{
    if (m.size() == _mi.size()) {
        for (size_t i = 0; i < m.size(); ++i) {
            _mi[i] = 1./m[i];
        }
    }
}
} // namespace nfm
