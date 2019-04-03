#include "nfm/ConjGrad.hpp"
#include "nfm/DynamicDescent.hpp"
#include "nfm/Adam.hpp"
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
    cout << "We will always start at (0, 3)." << endl << endl << endl;

    cout << "Let's first use CG to minimize the noiseless function." << endl;

    RosenbrockFunction<2> rbfun;
    double initpos[2] {0., 3.};

    LogManager::setLoggingFilePath("cg.out");
    ConjGrad cg(&rbfun);
    cg.setStepSize(0.01);
    cg.setX(initpos);
    cg.findMin();
    reportMinimum(cg);

    cout << "It took a bit longer than for the Parabola, but the minimum" << endl;
    cout << "was found to more than 10^-5 precision within 20 gradient steps." << endl << endl;
    cout << endl << "Now let's try some of the SGD optimizers." << endl << endl;
    cout << "We start with simple momentum SGD:" << endl;

    LogManager::setLoggingFilePath("sgdm.out");
    DynamicDescent dd(&rbfun);
    dd.setEpsX(0.); // don't stop on small x changes
    dd.setStepSize(0.001); // we need this small size to not overshoot too extremely
    dd.setMaxNIterations(1000);
    dd.setX(initpos);
    dd.findMin();
    reportMinimum(dd);

    cout << "That was the result after 1000(!) gradient steps," << endl;
    cout << "so simple momentum SGD doesn't work very well here." << endl << endl;
    cout << "Nesterov momentum:" << endl;

    LogManager::setLoggingFilePath("nest.out");
    dd.useNesterov();
    dd.setStepSize(0.001);
    dd.setBeta(0.95);
    dd.setX(initpos);
    dd.findMin();
    reportMinimum(dd);

    cout << "Nesterov does very good, but still takes more steps than CG." << endl << endl;

    cout << "AdaGrad:" << endl;

    LogManager::setLoggingFilePath("adag.out");
    dd.useAdaGrad();
    dd.setStepSize(1.31);
    dd.setX(initpos);
    dd.findMin();
    reportMinimum(dd);

    cout << "RMSProp:" << endl;

    LogManager::setLoggingFilePath("rmsp.out");
    dd.useRMSProp();
    dd.setStepSize(0.002); // this is a lucky stepsize
    dd.setBeta(0.95);
    dd.setX(initpos);
    dd.findMin();
    reportMinimum(dd);


    cout << "AdaDelta:" << endl;

    LogManager::setLoggingFilePath("adad.out");
    dd.useAdaDelta();
    dd.setStepSize(0.004);
    dd.setBeta(0.5);
    dd.setX(initpos);
    dd.findMin();
    reportMinimum(dd);

    cout << "Adam:" << endl;

    LogManager::setLoggingFilePath("adam.out");
    Adam adam(&rbfun, false, 0.09, 0.9, 0.9); // more parameters ftw
    adam.setEpsX(0.); // same as above
    adam.setMaxNIterations(1000);
    adam.setX(initpos);
    adam.findMin();
    reportMinimum(adam);


    cout << endl << "Now let's turn on the noise and try a combination of noisy CG and SGDM:" << endl;

    NoisyWrapper nrbf(&rbfun, 0.1); // sigma 0.1

    LogManager::setLoggingFilePath("cg-sgd_noise.out");
    ConjGrad cg2(&nrbf);

    // config
    cg2.setGradErrStop(false);
    cg2.setStepSize(0.2);
    cg2.setBackStep(0.1);
    cg2.useConjGradPR(); // when not stopping on rejection, better use Polak-Ribiere

    // stopping criteria
    cg2.setEpsX(0.);
    cg2.setEpsF(0.);
    cg2.setMaxNIterations(25); // at this point the method usually gets stuck due to noise

    cg2.setX(initpos);
    cg2.findMin();
    cout << "Init CG result:" << endl;
    reportMinimum(cg2);

    DynamicDescent dd2(&nrbf, DDMode::SGDM, false);
    dd2.setStepSize(0.0025);
    dd2.setBeta(0.95);
    dd2.setMaxNConstValues(0);
    dd2.setMaxNIterations(275);
    dd2.setX(cg2.getX());
    dd2.findMin();
    cout << "Final SGD result:" << endl;
    reportMinimum(dd2);

    cout << "This way we obtain decent final minima with less steps (max 25+275)," << endl;
    cout << "despite the combination of difficult function and noise." << endl << endl;

    cout << "But the added noise does not really hurt some(!) of the SGD algorithms." << endl;
    cout << "For example, if we use Adam right from the start (300 steps), the result is similar:" << endl;
    LogManager::setLoggingFilePath("adam_noise.out");
    Adam adam2(&nrbf, true, 0.09, 0.9, 0.9);
    adam2.setMaxNConstValues(0); // let's skip this check
    adam2.setMaxNIterations(300); // and just run 300 adam steps
    adam2.setX(initpos);
    adam2.findMin();
    reportMinimum(adam2);

    // end
    return 0;
}

