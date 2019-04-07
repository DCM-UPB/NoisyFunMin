#ifndef NFM_FIRE_HPP
#define NFM_FIRE_HPP

#include "nfm/NoisyFunMin.hpp"
#include "nfm/MDIntegrators.hpp"

#include <algorithm>

namespace nfm
{

// FIRE: Fast Inertial Relaxation Engine (https://doi.org/10.1103/PhysRevLett.97.170201)
// Implemented and extended by Jan Kessler
//
// This algorithm evolves a molecular dynamics (-like) trajectory of the model
// parameters in the potential energy given by the target function. It has a
// variable time step mechanism and will zero the velocity when going uphill.
//
// NOTE: By default, the original default FIRE algorithm is reproduced, using
//       Velocity-Verlet as MD integrator. But as extensions to the original,
//       we also provide the possibility to set per-index masses which will also
//       be acknowledged in the mixing term, you may use a more selective freezing
//       mechanism (freezing only the "offending" velocity elements) and finally you
//       have the ability to set a minimal time step together with a number of steps
//       to stay at that limit before terminating.
class FIRE: public NFM
{
protected:
    // the time-step parameters
    double _dtmax; // maximal MD time step
    double _dt0; // initial MD time step (default 0.1 * dtmax)
    double _dtmin = 0.; // minimal MD time step
    int _Ndtmin = 0; // number of steps with minimal dt until termination (disabled by default)

    // other FIRE parameters
    int _Nwait = 5; // number of steps with non-negative P until dt increase / alpha decrease
    double _finc = 1.1; // dt increase factor, el (1,inf)
    double _fdec = 0.5; // dt decrease factor, el (0,1)
    double _alpha0 = 0.1; // initial/reset mixing factor, el (0,1)
    double _falpha = 0.99; // alpha decrease factor, el (0,1)

    // MD / extension parameters
    md::Integrator _mdi = md::Integrator::VerletV; // MD integrator to be used, default Verlocity Verlet
    bool _flag_fullFreeze = true; // set to false if you want to use selective instead of global (original) freezing
    std::vector<double> _mi; // inverse masses (default all 1)

    // --- Internal methods
    bool _initializeMD(md::MDView &view, double dt);
    void _updateTarget(md::MDView &view, double dt);
    bool _isNDtMinReached(int Nmin);
    void _findMin() override;

public:
    explicit FIRE(NoisyFunctionWithGradient * targetfun, double dtmax, double dt0 = 0. /*will be set to 0.1*dtmax*/);
    ~FIRE() override = default;

    // Getters
    double getDt0() const { return _dt0; }
    double getDtMax() const { return _dtmax; }
    double getDtMin() const { return _dtmin; }
    int getNDtMin() const { return _Ndtmin; }

    int getNWait() const { return _Nwait; }
    double getFInc() const { return _finc; }
    double getFDec() const { return _fdec; }
    double getAlpha0() const { return _alpha0; }
    double getFAlpha() const { return _falpha; }

    md::Integrator getMDIntegrator() const { return _mdi; }
    bool getFullFreeze() const { return _flag_fullFreeze; }

    // Setters
    void setDt0(double dt0) { _dt0 = std::max(_dtmin, std::min(_dtmax, dt0)); }
    void setDtMax(double dtmax) { _dtmax = std::max(_dt0, dtmax); }
    void setDtMin(double dtmin) { _dtmin = std::max(0., std::min(_dt0, dtmin)); }
    void setNDtMin(int Ndtmin) { _Ndtmin = std::max(0, Ndtmin); }

    void setNWait(int Nwait) { _Nwait = std::max(0, Nwait); }
    void setFInc(double finc) { _finc = std::max(1., finc); }
    void setFDec(double fdec) { _fdec = std::max(0., std::min(1., fdec)); }
    void setAlpha0(double alpha0) { _alpha0 = std::max(0., std::min(1., alpha0)); }
    void setFAlpha(double falpha) { _falpha = std::max(0., std::min(1., falpha)); }

    void setMDIntegrator(md::Integrator mdi) { _mdi = mdi; }
    void setFullFreeze() { _flag_fullFreeze = true; }
    void setSelectiveFreeze() { _flag_fullFreeze = false; }
    void setMasses(const std::vector<double> &m);
    void resetMasses() { std::fill(_mi.begin(), _mi.end(), 1.); }
};
} // namespace nfm

#endif
