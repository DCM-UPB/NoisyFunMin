#ifndef ONED_TOOLS
#define ONED_TOOLS

#include "nfm/NoisyFunctionValue.hpp"
#include "nfm/Definitions.hpp"
#include "nfm/FunProjection1D.hpp"

#include <string>



namespace nfm
{

    void findBracket(NoisyFunction * f1d, NoisyFunctionValue &a, NoisyFunctionValue &b, NoisyFunctionValue &c);
    //               ^function   ^three points that will provide the bracket. a contains the starting point (i.e. is also an input)

    void parabgoldMinimization(NoisyFunction * f1d, const double &eps, NoisyFunctionValue &a, NoisyFunctionValue &b, NoisyFunctionValue &c);
    //                                  ^function   ^level of precision ^3 points that initially provide the bracket, in the end the minimum will be in b



    // PRIVATE FUNCTIONS

    void _writeabcInLog(const std::string &key, NoisyFunctionValue &a, NoisyFunctionValue &b, NoisyFunctionValue &c);
    // function used internally for writing on the log

    void _abortFindBracket();  // function called when the findBracket() routine is taking too many computations
}

#endif
