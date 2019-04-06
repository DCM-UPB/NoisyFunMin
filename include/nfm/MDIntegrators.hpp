#ifndef NFM_MDINTEGRATORS_HPP
#define NFM_MDINTEGRATORS_HPP

#include "nfm/NoisyFunction.hpp"

#include <functional>
#include <algorithm>
#include <numeric>

namespace nfm
{
namespace md
{

enum class Integrator
{
    EulerE, /*explicit Euler*/
    VerletV /*Velocity Verlet*/
};

// Struct used to present a set of vectors
// from the optimizers to MD integrators.
struct MDView
{
    const std::vector<double> &mi; // inverse masses
    std::vector<double> &x; // position
    std::vector<double> &v; // velocity
    std::vector<double> &a; // acceleration (F*mi)
    NoisyGradient &F; // force (passed to NoisyFun)
};

// --- Functions

// Given an enum and data, perform the step
NoisyValue doMDStep(Integrator mdi, NoisyFunctionWithGradient &efun, MDView &view, double dt);

// helper to compute acceleration
void computeAcceleration(MDView &view);


// MDIntegrator step functions are of the form:
// std::function<NoisyValue(NoisyFunctionWithGradient &efun,          MDView &view,  double dt)>
//                  ^returns new energy value         ^energy/force function ^in/out MDView ^time step

// Euler
NoisyValue ExplicitEulerIntegrator(NoisyFunctionWithGradient &efun, MDView &view, double dt);

// Standard Velocity-Verlet, 4 step version
NoisyValue VelocityVerletIntegrator(NoisyFunctionWithGradient &efun, MDView &view, double dt);
} // namespace md
} // namespace nfm

#endif
