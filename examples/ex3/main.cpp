#include "nfm/ConjGrad.hpp"
#include "nfm/DynamicDescent.hpp"
#include "nfm/Adam.hpp"
#include "nfm/FIRE.hpp"
#include "nfm/IRENE.hpp"
#include "nfm/LogManager.hpp"

#include <iostream>
#include <memory>

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

void minimize(nfm::NFM &optimizer, nfm::NoisyFunctionWithGradient &tfun, const std::vector<double> &initpos, const std::string &fname)
{
    nfm::LogManager::setLoggingFilePath(fname);
    optimizer.findMin(tfun, initpos);
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
    cout << "We will always start at (-0.1, 2.9).    " << endl;
    cout << "This ensures a formidable challenge!" << endl << endl;

    cout << "Let's first use Fletcher-Reeves CG to minimize the noiseless function:" << endl;

    RosenbrockFunction<2> rbfun;
    const std::vector<double> initpos{-0.1, 2.9};

    ConjGrad cg(rbfun.getNDim());
    cg.setStepSize(0.01);
    cg.disableStopping();
    cg.setMaxNIterations(30);
    minimize(cg, rbfun, initpos, "cgfr.out");

    cout << "It took a bit longer than for the Parabola, but the minimum" << endl;
    cout << "was found to about 10^-5 precision within 30 gradient steps." << endl << endl;
    cout << "Results for other SD and CG variants (after 30 gradient steps)..." << endl;

    cout << "Raw gradient/Steepest Descent:" << endl;
    cg.useRawGrad();
    minimize(cg, rbfun, initpos, "cgsd.out");

    cout << "Polak-Ribiere CG:" << endl;
    cg.useConjGradPR();
    minimize(cg, rbfun, initpos, "cgpr.out");

    cout << "Polak-Ribiere CG with reset:" << endl;
    cg.useConjGradPR0();
    minimize(cg, rbfun, initpos, "cgpr0.out");


    cout << endl << "Now let's try some of the SGD optimizers." << endl << endl;
    cout << "We start with simple momentum SGD:" << endl;

    DynamicDescent dd(rbfun.getNDim()); // default is SGD with momentum
    dd.setStepSize(0.0005); // despite the small step we still overshoot (but we move fast!)
    dd.setBeta(0.9);
    dd.disableStopping(); // we want only one stopping condition
    dd.setMaxNIterations(1000); // namely the fixed number of steps
    minimize(dd, rbfun, initpos, "sgdm.out");

    cout << "That was the result after 1000(!) gradient steps," << endl;
    cout << "so the SGD method is much slower than CG here (no noise!)." << endl << endl;
    cout << "Results for other SGD variants (1000 gradient steps):" << endl;
    cout << "Nesterov momentum:" << endl;

    dd.useNesterov();
    dd.setStepSize(0.001);
    dd.setBeta(0.9);
    minimize(dd, rbfun, initpos, "nest.out");


    cout << "RMSProp:" << endl;

    dd.useRMSProp();
    dd.setStepSize(0.004);
    dd.setBeta(0.9);
    minimize(dd, rbfun, initpos, "rmsp.out");


    cout << "AdaDelta:" << endl;

    dd.useAdaDelta();
    dd.setStepSize(0.004);
    dd.setBeta(0.9);
    minimize(dd, rbfun, initpos, "adad.out");

    cout << "Adam:" << endl;

    Adam adam(rbfun.getNDim(), false, 0.1);
    adam.setBeta1(0.9);
    adam.setBeta2(0.95); // more parameters ftw
    adam.disableStopping(); // same as above
    adam.setMaxNIterations(1000); // same as above
    minimize(adam, rbfun, initpos, "adam.out");


    cout << "And also the Fast Inertial Relaxation Engine (FIRE):" << endl;
    FIRE fire(rbfun.getNDim(), 0.05, 0.01);
    fire.disableStopping();
    fire.setMaxNIterations(1000);
    minimize(fire, rbfun, initpos, "fire.out");


    cout << endl << "Now let's turn on the noise and first try a combination of noisy CG and SGDM:" << endl;

    NoisyWrapper nrbf(&rbfun, 0.2); // sigma 0.2

    // config
    cg.setStepSize(0.02);
    cg.setBackStep(0.005);
    cg.setMaxNBracket(7);
    cg.useConjGradFR(); // in the noise case Fletcher-Reeves seems to do best

    // stopping criteria
    cg.disableStopping();
    cg.setMaxNIterations(30); // at this point the method usually gets stuck due to noise

    cout << "Init CG result:" << endl;
    minimize(cg, nrbf, initpos, "cg-sgd_noise.out");

    dd.useSGDM();
    dd.setStepSize(0.001);
    dd.setBeta(0.9);
    dd.disableStopping();
    dd.setMaxNIterations(470);

    cout << "Final SGD result:" << endl;
    minimize(dd, nrbf, cg.getX(), "cg-sgd_noise.out");

    cout << "This way we obtain decent final minima with less steps (max 30+470)," << endl;
    cout << "despite the combination of difficult function and noise." << endl << endl;

    cout << "But the added noise does not really hurt some(!) of the SGD algorithms." << endl;
    cout << "For example, if we use Adam right from the start (500 steps)," << endl;
    cout << "the result is very accurate if we turn averaging on:" << endl;

    adam.setAveraging(true);
    adam.setAlpha(0.1);
    adam.setBeta1(0.9);
    adam.setBeta2(0.95);
    adam.disableStopping();
    adam.setMaxNIterations(500); // just run 500 adam steps
    minimize(adam, nrbf, initpos, "adam_noise.out");

    cout << endl;
    cout << "Now with noise on and a shallow minimum, default FIRE got more problems:" << endl;
    FIRE fire2(nrbf.getNDim(), 0.03, 0.01);
    //fire2.setDtMin(0.01);
    //fire2.setSelectiveFreeze();
    //fire2.setNWait(0);
    fire2.disableStopping();
    fire2.setMaxNIterations(500);
    minimize(fire2, nrbf, initpos, "fire_noise.out");

    cout << "Our custom IRENE to the rescue ( same/default settings ):" << endl;
    IRENE irene(nrbf.getNDim(), 0.03, 0.01);
    //irene.setDtMin(0.005);
    //irene.setSelectiveFreeze();
    //irene.setNWait(0);
    //irene.setBeta(0.9);
    irene.disableStopping();
    irene.setMaxNIterations(500);
    minimize(irene, nrbf, initpos, "irene_noise.out");

    cout << "And with slightly optimized settings:" << endl;
    NoisyValue::setSigmaLevel(1.);
    IRENE irene2(nrbf.getNDim(), 0.03, 0.01);
    //irene2.setDtMin(0.001);
    //irene2.setSelectiveFreeze();
    //irene2.setNWait(0);
    irene2.setBeta(0.9);
    irene2.disableStopping();
    irene2.setMaxNIterations(500);
    minimize(irene2, nrbf, initpos, "irene2_noise.out");

    //cout << "And with beta>0:" << endl;
    //sirene.setBeta(0.5);
    //minimize(irene, initpos, "irene-b_noise.out");

    // end
    return 0;
}

