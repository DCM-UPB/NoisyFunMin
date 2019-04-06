#ifndef NFM_SIRENE_HPP
#define NFM_SIRENE_HPP

#include "nfm/FIRE.hpp"

namespace nfm
{

// Smooth Inertial Relaxation Engine for Noisy Energy surfaces (SIRENE)
// Author: Jan Kessler (2019)
//
// Inspired by FIRE (see FIRE.hpp or https://doi.org/10.1103/PhysRevLett.97.170201),
// this modified algorithm suffers less from the presence of significant statistical
// noise in the target function (energy) and corresponding gradients.
class SIRENE: public FIRE // reuse some members from FIRE
{
private:
    double _beta = 0.; // builds up averaged gradient to use for scalar product P (0 means use pure raw grad)
    bool _flag_soft = false; // don't stop completely on negative scalar product

    // --- Internal methods
    void _findMin() override;

public:
    explicit SIRENE(NoisyFunctionWithGradient * targetfun, double dtmax, double dt0): FIRE(targetfun, dtmax, dt0) {}

    // Getters
    double getBeta() const { return _beta; }
    bool getSoftFreeze() const { return _flag_soft; }

    // Setters
    void setBeta(double beta) { _beta = std::max(0., std::min(1., beta)); }
    void setSoftFreeze() { _flag_soft = true; }
    void setHardFreeze() { _flag_soft = false; }
};
} // namespace nfm

#endif
