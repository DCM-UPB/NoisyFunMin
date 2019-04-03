#include "nfm/ConjGrad.hpp"
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

    TestParabola3D nlp;
    std::vector<double> initpos{2.5, 1., -1.};
    auto cg = std::make_unique<ConjGrad>(&nlp);

    cg->setX(initpos);
    cg->findMin();
    reportMinimum(*cg);

    cout << "Now we do the same, but with noise added to the function and its gradient." << endl << endl;

    cout << "First we try a small amount of noise:" << endl;
    NoisyWrapper np1(&nlp, 0.05); // sigma = 0.05
    cg = std::make_unique<ConjGrad>(&np1);

    cg->setX(initpos);
    cg->findMin();
    reportMinimum(*cg);

    cout << "Now we try a larger amount of noise:" << endl;
    NoisyWrapper np2(&nlp, 0.25); // sigma = 0.25
    cg = std::make_unique<ConjGrad>(&np2);

    // Let's change some settings for this very noisy CG
    cg->setStepSize(2.);
    cg->setBackStep(1.); // allow small backstep for bracketing

    cg->setX(initpos);
    cg->findMin();
    reportMinimum(*cg);

    cout << endl << "As we can see, the CG method usually still converges to the true minimum,";
    cout << endl << "within the error due to noise. But the method will not improve the";
    cout << endl << "result any further. Hence, in many cases regular stochastic optimizers";
    cout << endl << "will achieve a better final result (but not in just 2 steps!)." << endl;
    cout << endl << "But if you have the possibility to decrease target function noise on demand";
    cout << endl << "you may try something like the following example using a custom policy.";
    cout << endl << "Note that this is just an artifical demonstration of the technique!" << endl << endl;

    // let's define a policy lambda (not binding the target function above directly)
    const double scale = 1./sqrt(2.); // supposing we double the number of samples each iteration, so error would go down like this
    auto myPolicy = [=](NFM &nfm, NoisyFunction &targetfun) {
        // decrease sigma iteratively
        auto * myfun = dynamic_cast<NoisyWrapper *>(&targetfun); // a static cast would be faster, but this is safer
        if (myfun) { myfun->setSigma(scale*myfun->getSigma()); } // supposing that we have the option to decrease the noise (e.g. by sampling twice as much)

        // we may also change settings of NFM/CG here
        auto * mycg = dynamic_cast<ConjGrad *>(&nfm);
        if (mycg) {
            // this checks if the step was rejected (usually due to impossible noisy bracketing)
            if (mycg->getDeltaX() <= 0.) { // then increase initial bracket
                mycg->setStepSize(2.*mycg->getStepSize());
                mycg->setBackStep(2.*mycg->getBackStep());
            } // this is not stable for real use, because the size may only go up!
        }
    };

    // set your own policy
    cg->setPolicy(myPolicy);

    // and we change stopping criteria
    cg->setGradErrStop(false);
    cg->setEpsX(0.);
    cg->setEpsF(0.);
    cg->setMaxNIterations(10); // fixed amount of steps

    // let's try again
    cg->setX(initpos);
    cg->findMin();
    reportMinimum(*cg);

    cout << "NOTE: You may enable detailed logging by uncommenting" << endl;
    cout << "      one of two lines in the beginning of the example." << endl;

    // end
    return 0;
}
