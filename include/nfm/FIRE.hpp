#ifndef NFM_FIRE_HPP
#define NFM_FIRE_HPP

#include "nfm/NoisyFunMin.hpp"

namespace nfm
{

// Fast Inertial Relaxation Engine (FIRE)
// Paper: https://doi.org/10.1103/PhysRevLett.97.170201
// Or freely available version: https://www.math.uni-bielefeld.de/~gaehler/papers/fire.pdf
//
// This algorithm evolves a molecular dynamics (-like) trajectory of the model
// parameters in the potential energy given by the target function. It has a
// variable time step mechanism and will zero the velocity when going uphill.
//
class FIRE: public NFM
{
protected:
    double _dtmax; // maximal MD time step
    double _dt0; // initial MD time step (default 0.1 * dtmax)
    int _Nmin = 5; // number of steps with non-negative P until dt increase / alpha decrease
    double _finc = 1.1; // dt increase factor, el (1,inf)
    double _fdec = 0.5; // dt decrease factor, el (0,1)
    double _alpha0 = 0.1; // initial/reset mixing factor, el (0,1)
    double _falpha = 0.99; // alpha decrease factor, el (0,1)


    // --- Internal methods
    void _updateTarget();
    void _doMDStep(std::vector<double> &v, double dt);
    void _findMin() override;

public:
    explicit FIRE(NoisyFunctionWithGradient * targetfun, double dtmax, double dt0 = 0. /*will be set to 0.1*dtmax*/);
    ~FIRE() override = default;

    // Getters
    double getDt0() const { return _dt0; }
    double getDtMax() const { return _dtmax; }
    int getNMin() const { return _Nmin; }
    double getFInc() const { return _finc; }
    double getFDec() const { return _fdec; }
    double getAlpha0() const { return _alpha0; }
    double getFAlpha() const { return _falpha; }

    // Setters
    void setDt0(double dt0) { _dt0 = std::max(0., std::min(_dtmax, dt0)); }
    void setDtMax(double dtmax) { _dtmax = std::max(0., dtmax); }
    void setNMin(int Nmin ) { _Nmin = std::max(0, Nmin); }
    void setFInc(double finc) { _finc = std::max(1., finc); }
    void setFDec(double fdec) { _fdec = std::max(0., std::min(1., fdec)); }
    void setAlpha0(double alpha0) { _alpha0 = std::max(0., std::min(1., alpha0)); }
    void setFAlpha(double falpha) { _falpha = std::max(0., std::min(1., falpha)); }
};
} // namespace nfm

#endif
