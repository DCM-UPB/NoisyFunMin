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
//   - findBracket : Find a suitable bracket a,b,c, such that f(a)>f(b)
//                   and f(c)>(fb), while a<b<c.
//
//   - brentMinimization: Find x such that f(x) is minimal,
//                        given a valid initial bracket.
//
//   - multiLineMinimization: Uses FunProjection1D and algorithms above
//                            to minimize a multi-dimensional NoisyFunction
//                            along a line defined by p0 and dir.

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

void writeBracketToLog(const std::string &key, const NoisyBracket &bracket);
// function used for writing NoisyBracket to the log

// find and initial bracket
NoisyBracket findBracket(NoisyFunction &f1d, double initX1, double initX2);
//                                     ^1D function ^the starting points^

// brent minimization with noisy values
NoisyIOPair1D brentMinimization(NoisyFunction &f1d, NoisyBracket bracket, double eps = 1.e-8);
//                                            ^1D function       ^init bracket   ^level of precision

// helper to perform line-minimization of multi-dim function
NoisyIOPair multiLineMinimization(NoisyFunction &mdf, const std::vector<double> &p0, const std::vector<double> &dir, double eps, double initX1 = 0., double initX2 = 1.);
//                                              ^mult-dim fun                   ^n-dim point                  ^direction   ^precision   ^provide initial bracket^
} // namespace nfm


#endif
