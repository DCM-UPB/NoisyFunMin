#include <cassert>
#include <iostream>

#include "nfm/1DTools.hpp"
#include "nfm/ConjGrad.hpp"
#include "nfm/LogManager.hpp"

#include "TestNFMFunctions.hpp"


int main()
{
    using namespace std;
    using namespace nfm;

    LogManager log_manager;
    //log_manager.setLoggingOn();

    double f, df;

    // define 3 noisy function values
    NoisyValue p1(1);
    NoisyValue p2(1);
    NoisyValue p3(1);


    // check parabola   x^2   ...
    Parabola parabola;

    // ... starting from x=-1000
    p1.setX(-1000.);
    parabola.f(p1.getX(), f, df);
    p1.setF(f, df);
    nfm::findBracket(&parabola, p1, p2, p3);
    assert(p1.getX(0) < 0.);
    assert(p3.getX(0) > 0.);

    // ... starting from x=1000
    log_manager.logString("\n\n=========================================================================\n\n");
    p1.setX(1000.);
    parabola.f(p1.getX(), f, df);
    p1.setF(f, df);
    nfm::findBracket(&parabola, p1, p2, p3);
    assert(p1.getX(0) < 0.);
    assert(p3.getX(0) > 0.);

    // ... starting from x=0
    log_manager.logString("\n\n=========================================================================\n\n");
    p1.setX(0.);
    parabola.f(p1.getX(), f, df);
    p1.setF(f, df);
    nfm::findBracket(&parabola, p1, p2, p3);
    assert(p1.getX(0) < 0.);
    assert(p3.getX(0) > 0.);


    // check well function   -1 if (-1 < x < 1) else +1   ...
    Well well;

    // ... starting from x=-10
    log_manager.logString("\n\n=========================================================================\n\n");
    p1.setX(-10.);
    well.f(p1.getX(), f, df);
    p1.setF(f, df);
    nfm::findBracket(&well, p1, p2, p3);
    assert(p1.getX(0) < -1.);
    assert(p3.getX(0) > 1.);

    // ... starting from x=-1000
    log_manager.logString("\n\n=========================================================================\n\n");
    p1.setX(-1000.);
    well.f(p1.getX(), f, df);
    p1.setF(f, df);
    bool flag_exception_thrown = false;
    try {
        nfm::findBracket(&well, p1, p2, p3);
    }
    catch (exception &e) {
        flag_exception_thrown = true;
    }
    assert(flag_exception_thrown);

    // ... starting from x=10
    log_manager.logString("\n\n=========================================================================\n\n");
    p1.setX(10.);
    well.f(p1.getX(), f, df);
    p1.setF(f, df);
    nfm::findBracket(&well, p1, p2, p3);
    assert(p1.getX(0) < -1.);
    assert(p3.getX(0) > 1.);

    return 0;
}
