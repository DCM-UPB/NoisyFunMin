#include <cassert>
#include <iostream>

#include "nfm/LineSearch.hpp"
#include "nfm/ConjGrad.hpp"
#include "nfm/LogManager.hpp"

#include "TestNFMFunctions.hpp"


int main()
{
    using namespace std;
    using namespace nfm;

    LogManager::setLogLevel(LogLevel::VERBOSE);
    assert(LogManager::isLoggingOn());
    assert(LogManager::isVerbose());
    LogManager::setLogLevel(LogLevel::NORMAL);
    assert(LogManager::isLoggingOn());
    assert(!LogManager::isVerbose());
    LogManager::setLogLevel(LogLevel::OFF);
    assert(!LogManager::isLoggingOn());
    assert(!LogManager::isVerbose());

    //LogManager::setLoggingOn(true); // uncomment if you actually want printout

    // define bracket of 3 noisy function input/output pairs
    NoisyBracket bracket{};
    NoisyIOPair1D &p1 = bracket.a;
    NoisyIOPair1D &p2 = bracket.b;
    NoisyIOPair1D &p3 = bracket.c;

    std::vector<double> inp(1); // helper vector

    // check pwr4   x^4   ...
    PowerFour pwr4;

    // ... using a=-3.   b=-2.   c=5.
    p1.x = inp[0] = -3.;
    p1.f = pwr4.f(inp);
    p2.x = inp[0] = -2.;
    p2.f = pwr4.f(inp);
    p3.x = inp[0] = 5.;
    p3.f = pwr4.f(inp);

    p2 = nfm::brentMin(pwr4, bracket, 20, 1e-5, 1e-10);
    assert(p2.x < 0.1);
    assert(p2.x > -0.1);
    assert(p2.f.value < 0.00001);


    return 0;
}
