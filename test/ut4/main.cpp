#include <cassert>
#include <cmath>
#include <iostream>

#include "nfm/1DTools.hpp"
#include "nfm/DynamicDescent.hpp"
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


    // test DynamicDescent
    DynamicDescent dyndesc(&f3d, 0.03);
    x[0] = -2.;
    x[1] = 1.0;
    x[2] = 0.0;
    dyndesc.setX(x);
    dyndesc.findMin();

    assert(fabs(dyndesc.getX(0) - 1.0) < 0.1);
    assert(fabs(dyndesc.getX(1) + 1.5) < 0.1);
    assert(fabs(dyndesc.getX(2) - 0.5) < 0.1);


    return 0;
}
