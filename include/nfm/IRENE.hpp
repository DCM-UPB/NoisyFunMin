#ifndef NFM_IRENE_HPP
#define NFM_IRENE_HPP

#include "nfm/FIRE.hpp"

namespace nfm
{

// Inertial Relaxation Engine for Noisy Energy surfaces (IRENE)
// Author: Jan Kessler (2019)
//
// Inspired by FIRE (see FIRE.hpp or https://doi.org/10.1103/PhysRevLett.97.170201),
// this modified algorithm suffers less from the presence of significant statistical
// noise in the gradient of the target function (energy).
class IRENE: public FIRE // reuse some members from FIRE
{
private:
    //double _beta = 0.; // builds up averaged gradient to use for scalar product P (0 means use pure raw grad)

    // --- Internal methods
    void _findMin() override;

public:
    explicit IRENE(NoisyFunctionWithGradient * targetfun, double dtmax, double dt0 = 0.): FIRE(targetfun, dtmax, dt0) {}

    // Getters
    //double getBeta() const { return _beta; }

    // Setters
    //void setBeta(double beta) { _beta = std::max(0., std::min(1., beta)); }
};
} // namespace nfm

#endif
