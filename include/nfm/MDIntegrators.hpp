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
    std::vector<double> &x; // position
    std::vector<double> &v; // velocity
    std::vector<double> &a; // acceleration (F*mi)
    std::function<void()> &update; // force update callback function (usually a lambda)
};

// --- Functions

// MDIntegrator step functions are of the form:
// std::function<void(MDView &view,  double dt)>
//                           ^in/out MDView ^time step
//
// Starting from the previous step's information in MDView x,v,a ,
// they will perform a MD time step according to dt and call the
// provided update() callback to ask for calculation of updated
// acceleration values, to be stored in a.

// Euler
void ExplicitEulerIntegrator(MDView &view, double dt);

// Standard Velocity-Verlet, 4 step version
void VelocityVerletIntegrator(MDView &view, double dt);


// calls the right integrator according to enum
void doMDStep(Integrator mdi, MDView &view, double dt);


// helper to compute accelerations
inline void computeAcceleration(const std::vector<double> &F, const std::vector<double> &mi, std::vector<double> &a)
{
    std::transform(F.begin(), F.end(), mi.begin(), a.begin(), std::multiplies<>());
}
} // namespace md
} // namespace nfm

#endif
