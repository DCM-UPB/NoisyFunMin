#include <cassert>
#include <cmath>
#include <iostream>

#include "nfm/1DTools.hpp"
#include "nfm/DynamicDescent.hpp"

#include "TestNFMFunctions.hpp"


int main()
{
    using namespace std;
    using namespace nfm;

    //LogManager::setLoggingOn(true);

    // define 3D function that I want to minimise
    F3D f3d;
    // introduce array with the initial position
    double x[3]{-2., 1., 0.};

    // test DynamicDescent
    DynamicDescent dyndesc(&f3d);
    dyndesc.setEpsX(0.001);
    dyndesc.setEpsF(1.e-5);
    dyndesc.setX(x);
    dyndesc.findMin();

    assert(fabs(dyndesc.getX(0) - 1.0) < 0.1);
    assert(fabs(dyndesc.getX(1) + 1.5) < 0.1);
    assert(fabs(dyndesc.getX(2) - 0.5) < 0.1);

    return 0;
}
