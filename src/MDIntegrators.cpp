#include "nfm/MDIntegrators.hpp"

#include <algorithm>
#include <numeric>

namespace nfm
{
namespace md
{

void ExplicitEulerIntegrator(MDView &view, const double dt)
{
    for (size_t i = 0; i < view.x.size(); ++i) {
        view.x[i] += dt*view.v[i];
        view.v[i] += dt*view.a[i];
    }
    view.update();
}

// Standard Velocity-Verlet, 4 step version
void VelocityVerletIntegrator(MDView &view, const double dt)
{
    const size_t ndim = view.x.size();
    const double hdt = 0.5*dt;
    for (size_t i = 0; i < ndim; ++i) {
        view.v[i] += hdt*view.a[i];
        view.x[i] += dt*view.v[i];
    }
    view.update();
    for (size_t i = 0; i < ndim; ++i) {
        view.v[i] += hdt*view.a[i];
    }
}

// Given an enum and data, perform the step
void doMDStep(const Integrator mdi, MDView &view, const double dt)
{
    switch (mdi) {
    case Integrator::EulerE:
        ExplicitEulerIntegrator(view, dt);

    case Integrator::VerletV:
        VelocityVerletIntegrator(view, dt);
    }
}
} // namespace md
} // namespace nfm