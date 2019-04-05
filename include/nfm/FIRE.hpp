#ifndef NFM_FIRE_HPP
#define NFM_FIRE_HPP

#include "nfm/NoisyFunMin.hpp"
#include "nfm/NoisyFunction.hpp"
#include "nfm/NoisyValue.hpp"

#include <list>

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
// NOTE: Per default the original algorithm is replicated. However, the setting
//       flag softFreeze (def false) and the beta parameter (def 0) allow to use
//       our own modifications intended for the use with noisy target functions.
class FIRE: public NFM
{
private:
    double _dt0; // initial MD time step
    double _dtmax; // maximal MD time step
    int _Nmin = 5; // number of steps with non-negative P until dt increase / alpha decrease
    double _finc = 1.1; // dt increase factor, el (1,inf)
    double _fdec = 0.5; // dt decrease factor, el (0,1)
    double _alpha0 = 0.1; // initial/reset mixing factor, el (0,1)
    double _falpha = 0.99; // alpha decrease factor, el (0,1)
    double _beta = 0.; // builds up averaged gradient to use for scalar product P (0 means standard FIRE)
    bool _flag_soft = false; // don't stop completely on negative scalar product (false means standard FIRE)


    // --- Internal methods
    void _updateTarget();
    void _doMDStep(std::vector<double> &v, double dt);
    void _findMin() override;

public:
    explicit FIRE(NoisyFunctionWithGradient * targetfun, double dtinit, double dtmax);
    ~FIRE() override = default;

    // Getters
    double getDt0() const { return _dt0; }
    double getDtMax() const { return _dtmax; }
    int getNMin() const { return _Nmin; }
    double getFInc() const { return _finc; }
    double getFDec() const { return _fdec; }
    double getAlpha0() const { return _alpha0; }
    double getFAlpha() const { return _falpha; }
    double getBeta() const { return _beta; }
    bool getSoftFreeze() const { return _flag_soft; }

    // Setters
    void setDt0(double dt0) { _dt0 = std::max(0., std::min(_dtmax, dt0)); }
    void setDtMax(double dtmax) { _dtmax = std::max(0., dtmax); }
    void setNMin(int Nmin ) { _Nmin = std::max(0, Nmin); }
    void setFInc(double finc) { _finc = std::max(1., finc); }
    void setFDec(double fdec) { _fdec = std::max(0., std::min(1., fdec)); }
    void setAlpha0(double alpha0) { _alpha0 = std::max(0., std::min(1., alpha0)); }
    void setFAlpha(double falpha) { _falpha = std::max(0., std::min(1., falpha)); }
    void setBeta(double beta) { _beta = std::max(0., std::min(1., beta)); }
    void setSoftFreeze() { _flag_soft = true; }
    void setHardFreeze() { _flag_soft = false; }
};
} // namespace nfm

#endif
