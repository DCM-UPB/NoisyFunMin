#include <cassert>
#include <cmath>
#include <iostream>

#include "nfm/DynamicDescent.hpp"
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
    double x[3]{-2., 1., 0.};


    // test DynamicDescent
    DynamicDescent dyndesc(&f3d);
    dyndesc.setX(x);
    dyndesc.findMin();

    assert(fabs(dyndesc.getX(0) - 1.0) < 0.1);
    assert(fabs(dyndesc.getX(1) + 1.5) < 0.1);
    assert(fabs(dyndesc.getX(2) - 0.5) < 0.1);


    // also with averaging
    dyndesc.setAveraging(true);
    dyndesc.setMaxNConstValues(5);
    dyndesc.setX(x);
    dyndesc.findMin();

    assert(fabs(dyndesc.getX(0) - 1.0) < 0.1);
    assert(fabs(dyndesc.getX(1) + 1.5) < 0.1);
    assert(fabs(dyndesc.getX(2) - 0.5) < 0.1);


    // AdaGrad
    dyndesc.useAdaGrad();
    dyndesc.setStepSize(3.);
    dyndesc.setX(x);
    dyndesc.findMin();

    assert(fabs(dyndesc.getX(0) - 1.0) < 0.1);
    assert(fabs(dyndesc.getX(1) + 1.5) < 0.15);
    assert(fabs(dyndesc.getX(2) - 0.5) < 0.15);


    // AdaDelta
    dyndesc.useAdaDelta();
    dyndesc.setStepSize(0.05);
    const double oldBeta = dyndesc.getBeta();
    dyndesc.setBeta(0.5); // I had problems when leaving this at default
    dyndesc.setX(x);
    dyndesc.findMin();
    dyndesc.setBeta(oldBeta); // set back the original value

    assert(fabs(dyndesc.getX(0) - 1.0) < 0.1);
    assert(fabs(dyndesc.getX(1) + 1.5) < 0.1);
    assert(fabs(dyndesc.getX(2) - 0.5) < 0.1);


    // RMSProp
    dyndesc.useRMSProp();
    dyndesc.setStepSize(0.05);
    dyndesc.setX(x);
    dyndesc.findMin();

    assert(fabs(dyndesc.getX(0) - 1.0) < 0.1);
    assert(fabs(dyndesc.getX(1) + 1.5) < 0.1);
    assert(fabs(dyndesc.getX(2) - 0.5) < 0.1);


    // Nesterov
    dyndesc.useNesterov();
    dyndesc.setStepSize(0.025);
    dyndesc.setX(x);
    dyndesc.findMin();

    assert(fabs(dyndesc.getX(0) - 1.0) < 0.1);
    assert(fabs(dyndesc.getX(1) + 1.5) < 0.1);
    assert(fabs(dyndesc.getX(2) - 0.5) < 0.1);

    return 0;
}
