#include <cassert>
#include <iostream>

#include "nfm/1DTools.hpp"
#include "nfm/ConjGrad.hpp"
#include "nfm/LogManager.hpp"

#include "TestNFMFunctions.hpp"

nfm::NoisyBracket prepareBracket(nfm::NoisyFunction &fun, const double ax, const double cx)
{
    nfm::NoisyBracket bracket{{ax,            {}},
                              {0.5*(ax + cx), {}},
                              {cx,            {}}}; // bracket from ax to cx
    std::vector<double> xvec(1);
    xvec[0] = bracket.a.x;
    bracket.a.f = fun(xvec);
    // we don't need to calculate the b function value
    xvec[0] = bracket.c.x;
    bracket.c.f = fun(xvec);
    return bracket;
}

int main()
{
    using namespace std;
    using namespace nfm;

    // test some log manager stuff
    LogManager::setLoggingOn(true);
    assert(LogManager::isLoggingOn());
    assert(LogManager::isVerbose());
    LogManager::setLoggingOn(false);
    assert(LogManager::isLoggingOn());
    assert(!LogManager::isVerbose());
    LogManager::setLoggingOff();
    assert(!LogManager::isLoggingOn());
    assert(!LogManager::isVerbose());

    //LogManager::setLoggingOn(true); // uncomment if you actually want printout

    // a noisy function input/output pair and a bracket
    NoisyIOPair p(1);
    NoisyBracket bracket{};
    NoisyIOPair1D &a = bracket.a;
    NoisyIOPair1D &c = bracket.c;
    bool flag_success;

    // check parabola   x^2   ...
    Parabola parabola;

    // ... starting from interval [-1000, 1]
    bracket = prepareBracket(parabola, -1000., 1.);
    flag_success = nfm::findBracket(parabola, bracket);
    assert(flag_success);
    assert(a.x < 0.);
    assert(c.x > 0.);

    // ... starting from interval [-5, 1000]
    LogManager::logString("\n\n=========================================================================\n\n");
    bracket = prepareBracket(parabola, 1000., -5.); // should be sorted in findBracket
    flag_success = nfm::findBracket(parabola, bracket);
    assert(flag_success);
    assert(a.x < 0.);
    assert(c.x > 0.);

    // ... starting from interval [-1.5,10]
    LogManager::logString("\n\n=========================================================================\n\n");
    bracket = prepareBracket(parabola, -1.5, 10.);
    flag_success = nfm::findBracket(parabola, bracket);
    assert(flag_success);
    assert(a.x < 0.);
    assert(c.x > 0.);


    // check well function   -1 if (-1 < x < 1) else +1   ...
    Well well;

    // ... starting from interval [-1.25, 3.5]
    LogManager::logString("\n\n=========================================================================\n\n");
    bracket = prepareBracket(well, -1.25, 5.5);
    flag_success = nfm::findBracket(well, bracket);
    assert(flag_success);
    assert(a.x < -1.);
    assert(c.x > 1.);

    // ... starting from interval [-1000, 500]
    LogManager::logString("\n\n=========================================================================\n\n");
    bracket = prepareBracket(well, -1000., 100.);
    flag_success = nfm::findBracket(well, bracket);
    assert(!flag_success);

    // ... starting from interval [-1, 3]
    LogManager::logString("\n\n=========================================================================\n\n");
    bracket = prepareBracket(well, -1.1, 3.);
    flag_success = nfm::findBracket(well, bracket);
    assert(flag_success);
    assert(a.x < -1.);
    assert(c.x > 1.);

    return 0;
}
