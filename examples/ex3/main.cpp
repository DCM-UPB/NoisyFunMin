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

    //LogManager::setLoggingOn(); // use this to enable log printout
    //LogManager::setLoggingOn(true); // use this for verbose printout of the method

    cout << endl;
    cout << "Testing Optimizers on Rosenbrock function" << endl << endl;
    cout << "We want to minimize the 2D function" << endl;
    cout << "    100*(y-x^2)^2 + (1-x)^2 ,      " << endl;
    cout << "which has its minimum in (1,1).    " << endl << endl;
    cout << "We will always start at (0, 3)." << endl << endl << endl;

    cout << "Let's first use CG to minimize the noiseless function." << endl;

    RosenbrockFunction<2> rbfun;
    double initpos[2] {0., 3.};

    ConjGrad cg(&rbfun);
    cg.setStepSize(0.01);
    cg.setX(initpos);
    cg.findMin();
    reportMinimum(cg);

    cout << "It took a bit longer than for the Parabola, but the minimum" << endl;
    cout << "was found to more than 10^-5 precision within 20 gradient steps." << endl << endl;
    cout << endl << "Now let's try some of the SGD optimizers." << endl << endl;
    cout << "We start with simple momentum SGD:" << endl;

    DynamicDescent dd(&rbfun);
    dd.setStepSize(0.001); // we need this small size to not overshoot initially
    dd.setMaxNIterations(1000);
    dd.setX(initpos);
    dd.findMin();
    reportMinimum(dd);

    cout << "That was the result after 1000(!) gradient steps," << endl;
    cout << "so simple momentum SGD doesn't work very well here." << endl << endl;
    cout << "Nesterov momentum:" << endl;

    dd.useNesterov();
    dd.setX(initpos);
    dd.findMin();
    reportMinimum(dd);

    cout << "It does better, but still takes 1000 steps to nearly match CG." << endl << endl;

    cout << "AdaDelta (poor):" << endl;

    dd.useAdaDelta();
    dd.setStepSize(0.004);
    dd.setBeta(0.5);
    dd.setX(initpos);
    dd.findMin();
    reportMinimum(dd);

    cout << "Adam (not good):" << endl;

    Adam adam(&rbfun, false, 0.008, 0.7, 0.9); // more parameters ftw
    adam.setMaxNIterations(1000);
    adam.setX(initpos);
    adam.findMin();
    reportMinimum(adam);


    cout << endl << "Now let's turn on the noise and try a combination of noisy CG and Adam:" << endl;

    NoisyWrapper nrbf(&rbfun, 0.05); // sigma 0.05
    ConjGrad cg2(&nrbf);

    // config
    cg2.setStepSize(0.2);
    cg2.setBackStep(0.1);
    cg2.useConjGradPR(); // when not stopping on rejection, better use Polak-Ribiere

    // stopping criteria
    cg2.setEpsX(0.);
    cg2.setEpsF(0.);
    cg2.setMaxNIterations(100);

    cg2.setX(initpos);
    cg2.findMin();
    cout << "Init CG result:" << endl;
    reportMinimum(cg2);

    Adam adam2(&nrbf, true, 0.1, 0.9, 0.9); // more parameters ftw
    adam2.setMaxNConstValues(0); // let's skip this check
    adam2.setMaxNIterations(200); // and just run 200 adam steps
    adam2.setX(cg2.getX()); // set result of CG
    adam2.findMin();
    cout << "Final Adam result:" << endl;
    reportMinimum(adam2);

    cout << "This way we obtain decent final minima with less steps (max 100+200)," << endl;
    cout << "despite the combination of difficult function and noise." << endl << endl << endl;

    // end
    return 0;
}

