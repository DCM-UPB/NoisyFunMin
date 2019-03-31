#include "nfm/ConjGrad.hpp"
#include "nfm/LogManager.hpp"

#include <iostream>

#include "../common/ExampleFunctions.hpp"

int main()
{
    using namespace std;
    using namespace nfm;

    cout << "We want to minimize the 3D function" << endl;
    cout << "    x^2 + (y+1)^2 + (z-2)^2" << endl;
    cout << "whose min is in (0, -1, 2)." << endl << endl;
    cout << "We will always start at (2.5, 1, -1)." << endl << endl;

    LogManager::setLoggingOn(); // use this to enable log printout
    //LogManager::setLoggingOn(true); // use this for verbose printout of the CG method

    cout << "We first minimize it, supposing to have no noise at all" << endl;

    Noiseless3DParabola nlp;
    auto cg = new ConjGrad(&nlp);

    std::vector<double> initpos{2.5, 1., -1.};
    cg->setX(initpos);

    cg->findMin();

    cout << "The found minimum is: ";
    cout << cg->getX(0) << "    " << cg->getX(1) << "    " << cg->getX(2) << endl << endl;
    delete cg;

    cout << "Now we repeat the minimisation adding a noise to the function and its gradient." << endl;

    cout << "First we try a small amount of noise:" << endl << endl;
    NoisyWrapper np1(&nlp, 0.05); // sigma = 0.05
    cg = new ConjGrad(&np1);
    cg->setX(initpos);

    cg->findMin();
    cout << "The found minimum is: ";
    cout << cg->getX(0) << "    " << cg->getX(1) << "    " << cg->getX(2) << endl << endl;
    delete cg;

    cout << "Now we try a larger amount of noise:" << endl << endl;
    NoisyWrapper np2(&nlp, 0.25); // sigma = 0.25
    cg = new ConjGrad(&np2);

    // Let's change some settings for this very noisy CG
    cg->setStepSize(0.5);
    cg->setBackStep(0.25); // allow small backstep for bracketing
    cg->setX(initpos);

    cg->findMin();
    cout << "The found minimum is: ";
    cout << cg->getX(0) << "    " << cg->getX(1) << "    " << cg->getX(2) << endl << endl;
    delete cg;

    cout << endl << "As we can see, the CG method still converges to the true minimum,";
    cout << endl << "within the error due to noise. But the method will not improve the";
    cout << endl << "result any further. Hence, in many cases regular stochastic optimizers";
    cout << endl << "will achieve a better final result (but not in just 2 steps!)." << endl;

    // end
    return 0;
}
