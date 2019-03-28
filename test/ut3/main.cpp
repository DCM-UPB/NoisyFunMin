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
    F3D f3d;
    // introduce array with the initial position
    double x[3];


    // test ConjGrad
    ConjGrad cjgrad(&f3d);
    x[0] = -2.;
    x[1] = 1.0;
    x[2] = 0.0;
    cjgrad.setX(x);
    cjgrad.setEpsF(1.e-5);
    cjgrad.findMin();
    assert(fabs(cjgrad.getX(0) - 1.0) < 0.15);
    assert(fabs(cjgrad.getX(1) + 1.5) < 0.15);
    assert(fabs(cjgrad.getX(2) - 0.5) < 0.1);

    return 0;
}
