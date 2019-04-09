#include <cassert>
#include <iostream>

#include "nfm/FIRE.hpp"
#include "nfm/IRENE.hpp"
#include "nfm/LogManager.hpp"

#include "TestNFMFunctions.hpp"


int main()
{
    using namespace std;
    using namespace nfm;

    LogManager::setLoggingOff();
    //LogManager::setLogLevel(LogLevel::VERBOSE);
    //LogManager::setLogLevel(LogLevel::NORMAL);

    // define 3D function that I want to minimise
    F3D f3d;
    std::vector<double> initpos{-2., 1., 0.};

    // --- Test default FIRE

    // with default integrator (Verlet)
    FIRE fire(&f3d, 1.);
    fire.setX(initpos);
    fire.findMin();

    assert(fabs(fire.getX(0) - 1.0) < 0.05);
    assert(fabs(fire.getX(1) + 1.5) < 0.05);
    assert(fabs(fire.getX(2) - 0.5) < 0.05);

    // with Euler integrator
    fire.setMDIntegrator(md::Integrator::EulerE);
    fire.setEpsX(0.); // euler doesn't move once when FIRE freezes system
    fire.setEpsF(0.);
    fire.setX(initpos);
    fire.findMin();

    assert(fabs(fire.getX(0) - 1.0) < 0.05);
    assert(fabs(fire.getX(1) + 1.5) < 0.05);
    assert(fabs(fire.getX(2) - 0.5) < 0.05);

    // --- Test extensions

    // test dtmin stopping
    fire.setMDIntegrator(md::Integrator::VerletV); // set back to original
    fire.disableStopping(); // we want to test the FIRE-specific stopping criterion
    fire.setDtMax(0.25);
    fire.setDt0(0.2);
    fire.setDtMin(0.15);
    fire.setNDtMin(1); // we stop when 1 step had dtmin time step
    fire.setX(initpos);
    fire.findMin();

    // which is way too early
    assert(fabs(fire.getX(0) - 1.0) > 0.05);
    assert(fabs(fire.getX(1) + 1.5) > 0.05);
    assert(fabs(fire.getX(2) - 0.5) > 0.05);

    // try selective freeze (and reenable fake const-list stopping)
    fire.setDtMax(1.);
    fire.setDtMin(0.);
    fire.setDt0(0.1);
    fire.setNDtMin(0);
    fire.setMaxNConstValues(20);
    fire.setSelectiveFreeze();
    fire.setX(initpos);
    fire.findMin();

    assert(fabs(fire.getX(0) - 1.0) < 0.05);
    assert(fabs(fire.getX(1) + 1.5) < 0.05);
    assert(fabs(fire.getX(2) - 0.5) < 0.05);

    // set different masses
    std::vector<double> m{0.75, 1.1, 0.9};
    fire.setMasses(m);
    fire.setX(initpos);
    fire.findMin();

    assert(fabs(fire.getX(0) - 1.0) < 0.05);
    assert(fabs(fire.getX(1) + 1.5) < 0.05);
    assert(fabs(fire.getX(2) - 0.5) < 0.05);

    // reset settings
    fire.setFullFreeze();
    fire.resetMasses();
    fire.setX(initpos);
    fire.findMin();

    // check that we properly reconfigured default state
    FIRE fire2(&f3d, 1.);
    fire2.setX(initpos);
    fire2.findMin();

    assert(fire.getX(0) == fire2.getX(0));
    assert(fire.getX(1) == fire2.getX(1));
    assert(fire.getX(2) == fire2.getX(2));

    // now we test also IRENE default (is identical here despite fake sigmas)
    IRENE irene(&f3d, 1.);
    irene.setX(initpos);
    irene.findMin();

    assert(irene.getX(0) == fire2.getX(0));
    assert(irene.getX(1) == fire2.getX(1));
    assert(irene.getX(2) == fire2.getX(2));
    assert(irene.getX() == fire2.getX());


    return 0;
}
