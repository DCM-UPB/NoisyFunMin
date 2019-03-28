#ifndef NFM_1DTOOLS_HPP
#define NFM_1DTOOLS_HPP

#include "nfm/FunProjection1D.hpp"
#include "nfm/NoisyValue.hpp"

#include <string>
#include <vector>

namespace nfm
{

// Here we provide functions to aid searching minima of NoisyFunctions
// along a line. If your NoisyFunction has multi-dimensional input,
// transform it to an effective 1D function by using FunProjection1D.
// The well-known line-search algorithms are implemented in typical form,
// but use the special overloads of NoisyValues.
//
// Functions:
//   - findBracket : Find a suitable bracket a,b,c, such that f(a)>f(b) and f(c)>(fb), while a<b<c.
//     Note 1: Adapted from GNU Scientific Libraries's bracketing code ( gsl/min/bracketing.c ),
//             mainly differing in making use of our NoisyValue comparison overloads.
//     Note 2: The method is not guaranteed to succeed, in pathologic cases even when there is
//             actually a minimum in the given initial interval. Check for the returned boolean.
//
//   - brentMinimization: Find x such that f(x) is minimal, given a valid initial bracket.
//     Note: Adapted from GNU Scientific Libraries's Brent Minimization code ( gsl/min/brent.c ),
//           mainly differing in making use of our NoisyValue comparison overloads and adapted
//           precision criterion.
//
//   - multiLineMinimization: Uses FunProjection1D and algorithms above to minimize a
//                            multi-dimensional NoisyFunction along a line defined by p0 and dir.
//

// 1D version of NoisyIOPair
struct NoisyIOPair1D
{
    double x;
    NoisyValue f;
};

struct NoisyBracket
// holds 3 bracketing positions and
// corresponding function values
{
    NoisyIOPair1D a;
    NoisyIOPair1D b;
    NoisyIOPair1D c;
};

// Function used for writing NoisyBracket to the log
void writeBracketToLog(const std::string &key, const NoisyBracket &bracket);

// Find an initial bracket for line minimizers
bool findBracket(NoisyFunction &f1d, NoisyBracket &bracket, double epsx = 1.e-8);
// ^did we have success        ^1D function  ^in/out bracket (a.x < b.x < c.x) ^bracket size tol

// Brent minimization with noisy values
NoisyIOPair1D brentMinimization(NoisyFunction &f1d, NoisyBracket bracket, double epsx = 1.e-8, double epsf = 1.e-8);
// ^minimized 1D-IO Pair                      ^1D function       ^init bracket   ^bracket size tol    ^target precision

// Helper to perform line-minimization of multi-dim function
// Returns the previous state if bracketing was not successful
NoisyIOPair multiLineMinimization(NoisyFunction &mdf, NoisyIOPair p0Pair, const std::vector<double> &dir, double initWidth = 1., double epsx = 1.e-8, double epsf = 1.e-8);
// ^minimized IO Pair                           ^mult-dim fun     ^last value with point            ^direction   ^init bracket size     ^x change tol        ^f change tol
} // namespace nfm


#endif
