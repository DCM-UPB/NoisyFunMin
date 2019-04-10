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
//
// This is achieved by the combination of two techniques: First of all, we make use of
// user-provided estimates of gradient errors, to make more safe decisions within the
// freezing / time tep mechanism. Furthermore, if the parameter "beta" is set > 0., we
// will accumulate the moving average gradient and use a mix of raw gradient and averaged
// gradient for the MD dynamics. The amount of mix-in is decided by the Signal-To-Noise
// ratio of gradient values. This allows to retain the very stiff and reactive dynamics
// of the original optimizer, but gains the ability to progress when gradients are noisy.
//
class IRENE: public FIRE // reuse some members from FIRE
{
private:
    double _beta = 0.; // exponential averaging factor for averaged acceleration

    // --- Internal methods
    void _findMin() override;

public:
    explicit IRENE(NoisyFunctionWithGradient * targetfun, double dtmax, double dt0 = 0.): FIRE(targetfun, dtmax, dt0) {}

    // Getters
    double getBeta() const { return _beta; }

    // Setters
    void setBeta(double beta) { _beta = std::max(0., std::min(1., beta)); }
};
} // namespace nfm

#endif
