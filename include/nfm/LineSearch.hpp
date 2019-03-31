#ifndef NFM_1DTOOLS_HPP
#define NFM_1DTOOLS_HPP

#include "nfm/FunProjection1D.hpp"
#include "nfm/NoisyValue.hpp"

#include <string>
#include <vector>

namespace nfm
{

// Here we provide functions to aid searching minima of NoisyFunctions along a line.
// If your NoisyFunction has multi-dimensional input, transform it to an effective 1D
// function by using FunProjection1D, before using the 1D methods. We use variants of
// well-known line-search algorithms, with some modifications and NoisyValue overloads.
//
// Functions:
//   - findBracket: Find a suitable bracket a,b,c, such that f(a)>f(b) and f(c)>(fb), while a<b<c.

//     Note 1: Originally based on GNU Scientific Libraries's bracketing code ( gsl/min/bracketing.c ),
//             but using NoisyValue overloads and other modifications.
//     Note 2: The algorithm requires an initial "bracket" [(a,fa),(b,fb),(c,fc)] with xa < xb < xc.
//             The bracket might be increased to the right side, but the initial lower bound will
//             never be lowered (but may go up).
//     Note 3: The method is not guaranteed to succeed, in the worst case even when there is
//             actually a minimum in the given initial interval. Check for the returned boolean.
//             If the boolean is true, the bracket is valid to be used for brentMin.
//
//
//   - brentMin: Find x such that f(x) is minimal, given a valid initial bracket.
//
//     Note: Adapted from GNU Scientific Libraries's Brent Minimization code ( gsl/min/brent.c ),
//           differing in making use of our NoisyValue comparison overloads and other improvements
//           for the noisy minimization use case.
//
//
//   - multiLineMin: Uses FunProjection1D and findBracket/Brent to minimize a multi-dimensional
//                   NoisyFunction along a line defined by last point p0 and direction dir. Given
//                   left and right steps define the (initial) search interval around p0.
//
//     Note 1: As stated above, the left boundary is hard while the right one may be increased.
//     Note 2: We set a third initial point at the golden section point between left and right.
//             But if the left step is passed as 0, the known function value passed via p0Pair
//             will be used as function value for the lower boundary at 0 (to save an evaluation).
//

// Global values to use in 1D Algos
namespace m1d_detail
{
static constexpr double GOLDEN = 1.618033988749895; // Golden Ratio
static constexpr double IGOLD2 = 1/(GOLDEN*GOLDEN); // 0.38196601125010515
static constexpr double STD_XTOL = 1.e-5; // default x tolerance
static constexpr double STD_FTOL = 1.e-8; // default f tolerance
} // namespace m1d_detail


// --- Line-Search Structs

// 1D version of NoisyIOPair
struct NoisyIOPair1D
{
    double x;
    NoisyValue f;
};

// Holds 3 bracketing positions and function values
struct NoisyBracket
{
    NoisyIOPair1D a;
    NoisyIOPair1D b;
    NoisyIOPair1D c;
};

// Parameters for multi-dimensional line-search
struct MLMParams
{
    double stepLeft; // initial step left of p0 (may be 0)
    double stepRight; // initial step right of p0 (may not be 0)
    int maxNBracket; // maximal number of bracketing iterations
    int maxNMinimize; // maximal number of minimization iterations
    double epsx; // x distance tolerance
    double epsf; // f distance tolerance
};

inline MLMParams defaultMLMParams()
{
    return {.stepLeft = 0., .stepRight = 1.,
            .maxNBracket = 10, .maxNMinimize = 20,
            .epsx = m1d_detail::STD_XTOL, .epsf = m1d_detail::STD_FTOL};
}


// --- Line-Search Functions

// Function used for writing NoisyBracket to the log
void writeBracketToLog(const std::string &key, const NoisyBracket &bracket);

// Find a valid bracket (starts with bracket [A, B, C], may increase interval to the right)
bool findBracket(NoisyFunction &f1d, NoisyBracket &bracket, int maxNIter, double epsx = m1d_detail::STD_XTOL);
// ^did we have success        ^1D function  ^in/out bracket (a.x < b.x < c.x)   ^bracket size tol

// Brent minimization with noisy values, requires valid NoisyBracket with a.f > b.f, b.f < c.f and a.x < b.x < c.x
NoisyIOPair1D brentMin(NoisyFunction &f1d, NoisyBracket bracket, int maxNIter, double epsx = m1d_detail::STD_XTOL, double epsf = m1d_detail::STD_FTOL);
// ^minimized 1D-IO Pair             ^1D function       ^init bracket ^iter limit     ^bracket size tol                   ^target precision

// Helper to perform line-minimization of multi-dim function
// Returns the previous state if minimization was not successful
NoisyIOPair multiLineMin(NoisyFunction &mdf, NoisyIOPair p0Pair, const std::vector<double> &dir, MLMParams params = defaultMLMParams());
// ^minimized IO Pair                  ^multi-dim fun    ^last value with point            ^direction      ^configuration
} // namespace nfm

#endif
