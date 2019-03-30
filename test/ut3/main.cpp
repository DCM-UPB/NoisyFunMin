#include <cassert>
#include <cmath>
#include <iostream>

#include "nfm/ConjGrad.hpp"
#include "nfm/LogManager.hpp"

#include "TestNFMFunctions.hpp"


int main()
{
    using namespace std;
    using namespace nfm;

    LogManager::setLoggingOff();
    //LogManager::setLoggingOn(true);

    // define 3D function that I want to minimise
    F3D f3d;
    // introduce array with the initial position
    double x[3];
    x[0] = -2.;
    x[1] = 1.0;
    x[2] = 0.0;

    const double XTOL = 0.15;
    const double YTOL = 0.10;
    const double ZTOL = 0.10;


    // test ConjGrad (Fletcher-Reeves)
    ConjGrad cjgrad(&f3d);
    cjgrad.setX(x);
    cjgrad.findMin();
    assert(fabs(cjgrad.getX(0) - 1.0) < XTOL);
    assert(fabs(cjgrad.getX(1) + 1.5) < YTOL);
    assert(fabs(cjgrad.getX(2) - 0.5) < ZTOL);

    // steepest descent version
    cjgrad.useRawGrad();
    cjgrad.setX(x);
    cjgrad.findMin();
    assert(fabs(cjgrad.getX(0) - 1.0) < XTOL);
    assert(fabs(cjgrad.getX(1) + 1.5) < YTOL);
    assert(fabs(cjgrad.getX(2) - 0.5) < ZTOL);

    // Polak-Ribiere version
    cjgrad.useConjGradPR();
    cjgrad.setX(x);
    cjgrad.findMin();
    assert(fabs(cjgrad.getX(0) - 1.0) < XTOL);
    assert(fabs(cjgrad.getX(1) + 1.5) < YTOL);
    assert(fabs(cjgrad.getX(2) - 0.5) < ZTOL);

    // PR-CG with reset
    cjgrad.useConjGradPR0();
    cjgrad.setX(x);
    cjgrad.findMin();
    assert(fabs(cjgrad.getX(0) - 1.0) < XTOL);
    assert(fabs(cjgrad.getX(1) + 1.5) < YTOL);
    assert(fabs(cjgrad.getX(2) - 0.5) < ZTOL);

    return 0;
}
