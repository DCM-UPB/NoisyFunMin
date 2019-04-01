#include "nfm/ConjGrad.hpp"
#include "nfm/LogManager.hpp"

#include <iostream>

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
    cout << "Conjugate Gradient Example" << endl << endl;
    cout << "We want to minimize the 3D function" << endl;
    cout << "    x^2 + (y+1)^2 + (z-2)^2" << endl;
    cout << "whose min is in (0, -1, 2)." << endl << endl;
    cout << "We will always start at (2.5, 1, -1)." << endl << endl << endl;

    cout << "We first minimize the function, supposing to have no noise at all" << endl;

    Noiseless3DParabola nlp;
    std::vector<double> initpos{2.5, 1., -1.};
    auto cg = new ConjGrad(&nlp);

    cg->setX(initpos);
    cg->findMin();
    reportMinimum(*cg);
    delete cg;

    cout << "Now we do the same, but with noise added to the function and its gradient." << endl << endl;

    cout << "First we try a small amount of noise:" << endl;
    NoisyWrapper np1(&nlp, 0.05); // sigma = 0.05
    cg = new ConjGrad(&np1);

    cg->setX(initpos);
    cg->findMin();
    reportMinimum(*cg);
    delete cg;

    cout << "Now we try a larger amount of noise:" << endl;
    NoisyWrapper np2(&nlp, 0.25); // sigma = 0.25
    cg = new ConjGrad(&np2);

    // Let's change some settings for this very noisy CG
    cg->setStepSize(0.5);
    cg->setBackStep(0.25); // allow small backstep for bracketing

    cg->setX(initpos);
    cg->findMin();
    reportMinimum(*cg);
    delete cg;

    cout << endl << "As we can see, the CG method still converges to the true minimum,";
    cout << endl << "within the error due to noise. But the method will not improve the";
    cout << endl << "result any further. Hence, in many cases regular stochastic optimizers";
    cout << endl << "will achieve a better final result (but not in just 2 steps!)." << endl << endl;

    cout << "NOTE: You may enable detailed logging by uncommenting" << endl;
    cout << "      one of two lines in the beginning of the example." << endl;

    // end
    return 0;
}
