#include <cassert>
#include <cmath>
#include <iostream>

#include "nfm/1DTools.hpp"
#include "nfm/ConjGrad.hpp"
#include "nfm/LogManager.hpp"

#include "TestNFMFunctions.hpp"


int main()
{
    using namespace std;
    using namespace nfm;

    LogManager::setLoggingOn(true);

    // define 3D function that I want to minimise
    auto * f3d = new F3D();
    // introduce array with the initial position
    double x[3];


    // test ConjGrad
    ConjGrad cjgrad(f3d);
    x[0] = -2.;
    x[1] = 1.0;
    x[2] = 0.0;
    cjgrad.setX(x);
    cjgrad.findMin();
    assert(fabs(cjgrad.getX(0) - 1.0) < 0.1);
    assert(fabs(cjgrad.getX(1) + 1.5) < 0.1);
    assert(fabs(cjgrad.getX(2) - 0.5) < 0.1);

    delete f3d;

    return 0;
}
