#include <cassert>
#include <iostream>

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
    xvec[0] = bracket.b.x;
    bracket.b.f = fun(xvec);
    xvec[0] = bracket.c.x;
    bracket.c.f = fun(xvec);
    return bracket;
}

void assertBracket(const nfm::NoisyBracket &bracket, const double maxax, const double mincx)
{
    const nfm::NoisyIOPair1D &a = bracket.a;
    const nfm::NoisyIOPair1D &b = bracket.b;
    const nfm::NoisyIOPair1D &c = bracket.c;

    assert(a.x < maxax);
    assert(c.x > mincx);
    assert(a.x != b.x);
    assert(c.x != b.x);
    assert(a.f > b.f);
    assert(c.f > b.f);
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
    NoisyBracket bracket{};
    const int nbracket = 10; // max bracketing attempts
    bool flag_success;

    // check parabola   x^2   ...
    Parabola parabola;

    // ... starting from interval [-1000, -1]
    bracket = prepareBracket(parabola, -1000., -1.);
    flag_success = nfm::findBracket(parabola, bracket, nbracket);
    assert(flag_success);
    assertBracket(bracket, 0., 0.);

    // ... starting from interval [-5, 1000]
    LogManager::logString("\n\n=========================================================================\n\n");
    bracket = prepareBracket(parabola, 1000., -5.); // should be sorted in findBracket
    flag_success = nfm::findBracket(parabola, bracket, nbracket);
    assert(flag_success);
    assertBracket(bracket, 0., 0.);

    // ... starting from interval [-1.1,-0.9]
    LogManager::logString("\n\n=========================================================================\n\n");
    bracket = prepareBracket(parabola, -1.1, -0.9);
    flag_success = nfm::findBracket(parabola, bracket, nbracket);
    assert(flag_success);
    assertBracket(bracket, 0., 0.);


    // check well function   -1 if (-1 < x < 1) else +1   ...
    Well well;

    // ... starting from interval [-1.1, 0.9]
    LogManager::logString("\n\n=========================================================================\n\n");
    bracket = prepareBracket(well, -1.1, 0.9);
    flag_success = nfm::findBracket(well, bracket, nbracket);
    assert(flag_success);
    assertBracket(bracket, -1., 1.);

    // ... starting from interval [-1.1, -1]
    LogManager::logString("\n\n=========================================================================\n\n");
    bracket = prepareBracket(well, -1.1, -1.);
    flag_success = nfm::findBracket(well, bracket, nbracket);
    assert(flag_success);
    assertBracket(bracket, -1., 1.);

    // ... starting from interval [-5, -0]
    LogManager::logString("\n\n=========================================================================\n\n");
    bracket = prepareBracket(well, -5, 0.);
    flag_success = nfm::findBracket(well, bracket, nbracket);
    assert(flag_success);
    assertBracket(bracket, -1., 1.);

    // ... starting from interval [1, 2] ( this will fail )
    LogManager::logString("\n\n=========================================================================\n\n");
    bracket = prepareBracket(well, 1., 2.);
    flag_success = nfm::findBracket(well, bracket, nbracket);
    assert(!flag_success);

    // ... starting from interval [-1.5, 1.5] (should be fine immediately)
    LogManager::logString("\n\n=========================================================================\n\n");
    bracket = prepareBracket(well, -1.1, 1.1);
    flag_success = nfm::findBracket(well, bracket, nbracket);
    assert(flag_success);
    assertBracket(bracket, -1., 1.);

    return 0;
}
