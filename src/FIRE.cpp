#include "nfm/FIRE.hpp"

#include "nfm/LogManager.hpp"

#include <cmath>

namespace nfm
{

FIRE::FIRE(const int ndim, const double dtmax, const double dt0):
        NFM(ndim, true), _dtmax(std::max(0., dtmax)), _dt0((dt0 > 0.) ? std::min(dt0, dtmax) : 0.1*dtmax)
{
    _mi.assign(_grad.size(), 1.); // inverse masses default to 1
    // override defaults
    this->setEpsX(1.e-08); // FIRE moves can occasionally be very small (with Euler integration this check must be off)
    this->setGradErrStop(false); // don't stop on noisy-low gradients, by default
}

// --- Minimization

void FIRE::_findMin()
{
    LogManager::logString("\nBegin FIRE::findMin() procedure\n");

    // helper variables
    std::vector<double> v(_grad.size()); // velocity vector
    std::vector<double> a(_grad.size()); // acceleration vector (F*mi)
    std::function<void()> update = [&]()
    { // MD force update lambda
        _last.f = _gradfun->fgrad(_last.x, _grad);
        md::computeAcceleration(_grad.val, _mi, a);
    };
    md::MDView mdview{.x = _last.x, .v = v, .a = a, .update = update}; // references for MD integrator

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
            LogManager::logString("\nFIRE::findMin() Step " + std::to_string(iter) + "\n");
        }

        // compute the gradient and current target
        this->_updateTarget(mdview, dt);
        if (this->_isNDtMinReached(Nmin) || this->_shouldStop()) { break; } // we are done

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

    // compute initial step
    view.update(); // compute initial force and store acceleration
    for (int i = 0; i < _ndim; ++i) {
        view.v[i] += dt*view.a[i]; // we start with an initial velocity
    }
    // other stuff
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
    md::doMDStep(_mdi, view, dt);
    this->_storeLastValue();
    this->_writeGradientToLog();
}

bool FIRE::_isNDtMinReached(const int Nmin)
{
    if (_Ndtmin > 0 && Nmin > _Ndtmin) {
        LogManager::logString("\nStopping Reason: Maximal number of steps with minimal time step.\n");
        return true;
    }
    return false;
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
