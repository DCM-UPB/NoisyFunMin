#include "nfm/ConjGrad.hpp"
#include "nfm/DynamicDescent.hpp"
#include "nfm/Adam.hpp"
#include "nfm/FIRE.hpp"
#include "nfm/LogManager.hpp"

#include <iostream>
#include <memory>
#include <cassert>

#include "../common/ExampleFunctions.hpp"

void reportMinimum(const nfm::NFM &optimizer)
{
    using namespace std;
    cout << "The found minimum is: ";
    cout << "f(" << optimizer.getX(0);
    for (int i = 1; i < optimizer.getNDim(); ++i) {
        cout << ", " << optimizer.getX(i);
    }
    cout << ") = " << optimizer.getFDf() << endl << endl;
}

void minimize(nfm::NFM &optimizer, const std::vector<double> &initpos, const std::string &fname)
{
    nfm::LogManager::setLoggingFilePath(fname);
    optimizer.setX(initpos);
    optimizer.findMin();
    reportMinimum(optimizer);
}

int main()
{
    using namespace std;
    using namespace nfm;

    LogManager::setLoggingOn(true); // leave this on, because for plotting we are logging to files

    cout << endl;
    cout << "Testing Optimizers on Rosenbrock function" << endl << endl;
    cout << "We want to minimize the 2D function" << endl;
    cout << "    100*(y-x^2)^2 + (1-x)^2 ,      " << endl;
    cout << "which has its minimum in (1,1).    " << endl << endl;
    cout << "We will always start at (0, 3).    " << endl;
    cout << "This ensures a formidable challenge!" << endl << endl;

    cout << "Let's first use Fletcher-Reeves CG to minimize the noiseless function:" << endl;

    RosenbrockFunction<2> rbfun;
    const std::vector<double> initpos{0., 3.};

    ConjGrad cg(&rbfun);
    cg.setStepSize(0.01);
    cg.disableStopping();
    cg.setMaxNIterations(25);
    minimize(cg, initpos, "cgfr.out");

    cout << "It took a bit longer than for the Parabola, but the minimum" << endl;
    cout << "was found to more than 10^-8 precision within 25 gradient steps." << endl << endl;
    cout << "Results for other SD and CG variants..." << endl;

    cout << "Raw gradient/Steepest Descent:" << endl;
    cg.useRawGrad();
    minimize(cg, initpos, "cgsd.out");

    cout << "Polak-Ribiere CG:" << endl;
    cg.useConjGradPR();
    minimize(cg, initpos, "cgpr.out");

    cout << "Polak-Ribiere CG with reset:" << endl;
    cg.useConjGradPR0();
    minimize(cg, initpos, "cgpr0.out");


    cout << endl << "Now let's try some of the SGD optimizers." << endl << endl;
    cout << "We start with simple momentum SGD:" << endl;

    DynamicDescent dd(&rbfun); // default is SGD with momentum
    dd.setStepSize(0.001); // despite the small step we still overshoot (but we move fast!)
    dd.disableStopping(); // we want only one stopping condition
    dd.setMaxNIterations(1000); // namely the fixed number of steps
    minimize(dd, initpos, "sgdm.out");

    cout << "That was the result after 1000(!) gradient steps," << endl;
    cout << "so simple momentum SGD doesn't work very well here." << endl << endl;
    cout << "Nesterov momentum:" << endl;

    dd.useNesterov();
    dd.setStepSize(0.001);
    dd.setBeta(0.95);
    minimize(dd, initpos, "nest.out");

    cout << "Nesterov does very good, but still takes more steps than CG." << endl << endl;

    cout << "AdaGrad:" << endl;

    dd.useAdaGrad();
    dd.setStepSize(1.31);
    minimize(dd, initpos, "adag.out");

    cout << "RMSProp:" << endl;

    dd.useRMSProp();
    dd.setStepSize(0.002); // this is a lucky stepsize
    dd.setBeta(0.95);
    minimize(dd, initpos, "rmsp.out");


    cout << "AdaDelta:" << endl;

    dd.useAdaDelta();
    dd.setStepSize(0.004);
    dd.setBeta(0.5);
    minimize(dd, initpos, "adad.out");

    cout << "Adam:" << endl;

    Adam adam(&rbfun, false, 0.09);
    adam.setBeta1(0.9);
    adam.setBeta2(0.9); // more parameters ftw
    adam.disableStopping(); // same as above
    adam.setMaxNIterations(1000); // same as above
    minimize(adam, initpos, "adam.out");


    cout << endl << "Now let's turn on the noise and try a combination of noisy CG and SGDM:" << endl;

    NoisyWrapper nrbf(&rbfun, 0.1); // sigma 0.1
    ConjGrad cg2(&nrbf);

    // config
    cg2.setStepSize(0.2);
    cg2.setBackStep(0.1);
    cg2.useConjGradPR(); // when not stopping on rejection, better use Polak-Ribiere

    // stopping criteria
    cg2.disableStopping();
    cg2.setMaxNIterations(25); // at this point the method usually gets stuck due to noise

    cout << "Init CG result:" << endl;
    minimize(cg2, initpos, "cg-sgd_noise.out");

    DynamicDescent dd2(&nrbf, DDMode::SGDM, false);
    dd2.setStepSize(0.0025);
    dd2.setBeta(0.95);
    dd2.disableStopping();
    dd2.setMaxNIterations(275);

    cout << "Final SGD result:" << endl;
    minimize(dd2, cg2.getX(), "cg-sgd_noise.out");

    cout << "This way we obtain decent final minima with less steps (max 25+275)," << endl;
    cout << "despite the combination of difficult function and noise." << endl << endl;

    cout << "But the added noise does not really hurt some(!) of the SGD algorithms." << endl;
    cout << "For example, if we use Adam right from the start (300 steps), the result is similar:" << endl;

    Adam adam2(&nrbf, true, 0.09);
    adam2.setBeta1(0.9);
    adam2.setBeta2(0.9);
    adam2.disableStopping();
    adam2.setMaxNIterations(300); // just run 300 adam steps
    minimize(adam2, initpos, "adam_noise.out");

    cout << endl;
    cout << "Last but not least, here are the results for the (Soft) Fast Inertial Relaxation Engine (FIRE):" << endl;
    FIRE fire(&nrbf, 0.005, 0.02);
    //fire.setNMin(3);
    fire.setBeta(0.7);
    fire.setSoftFreeze();
    fire.disableStopping();
    fire.setMaxNIterations(300);
    minimize(fire, initpos, "fire_noise.out");

    // end
    return 0;
}

