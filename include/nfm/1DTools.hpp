#ifndef NFM_1DTOOLS_HPP
#define NFM_1DTOOLS_HPP

#include "nfm/FunProjection1D.hpp"
#include "nfm/NoisyValue.hpp"

#include <string>
#include <tuple>

namespace nfm
{

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

NoisyBracket findBracket(NoisyFunction &f1d, double initX);
//                                   ^1D function     ^the starting point

NoisyIOPair1D parabgoldMinimization(NoisyFunction &f1d, const double &eps, NoisyBracket bracket);
//                                               ^1D function     ^level of precision   ^provide initial bracket
} // namespace nfm

#endif
