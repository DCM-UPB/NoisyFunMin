#include "nfm/MDIntegrators.hpp"

#include <algorithm>
#include <numeric>

namespace nfm
{
namespace md
{

void computeAcceleration(MDView &view)
{
    std::transform(view.F.val.begin(), view.F.val.end(), view.mi.begin(), view.a.begin(), std::multiplies<>());
}

NoisyValue ExplicitEulerIntegrator(NoisyFunctionWithGradient &efun, MDView &view, const double dt)
{
    for (int i = 0; i < efun.getNDim(); ++i) {
        view.v[i] += dt*view.a[i];
        view.x[i] += dt*view.v[i];
    }
    const NoisyValue newE = efun.fgrad(view.x, view.F);
    computeAcceleration(view);
    return newE;
}

// Standard Velocity-Verlet, 4 step version
NoisyValue VelocityVerletIntegrator(NoisyFunctionWithGradient &efun, MDView &view, const double dt)
{
    const int ndim = efun.getNDim();
    for (int i = 0; i < ndim; ++i) {
        view.v[i] += 0.5*dt*view.a[i];
        view.x[i] += dt*view.v[i];
    }
    const NoisyValue newE = efun.fgrad(view.x, view.F);
    computeAcceleration(view);
    for (int i = 0; i < ndim; ++i) {
        view.v[i] += 0.5*dt*view.a[i];
    }
    return newE;
}

// Given an enum and data, perform the step
NoisyValue doMDStep(Integrator mdi, NoisyFunctionWithGradient &efun, MDView &view, const double dt)
{
    switch (mdi) {
    case Integrator::EulerE:
        return ExplicitEulerIntegrator(efun, view, dt);

    case Integrator::VerletV:
        return VelocityVerletIntegrator(efun, view, dt);

    default:
        return {0., 0.}; // cannot happen
    }
}
} // namespace md
} // namespace nfm