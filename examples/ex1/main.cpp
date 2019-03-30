#include "nfm/ConjGrad.hpp"
#include "nfm/LogManager.hpp"

#include <iostream>

#include "../common/ExampleFunctions.hpp"

int main()
{
    using namespace std;
    using namespace nfm;

    cout << "We want to minimize the 2D function" << endl;
    cout << "    (x-1)^2 + (y+2)^2" << endl;
    cout << "whose min is in (1, -2)." << endl << endl << endl;

    LogManager::setLoggingOn(); // use this to enable log printout
    //LogManager::setLoggingOn(true); // use this for verbose printout of the CG method

    cout << "we first minimize it, supposing to have no noise at all" << endl;

    Noiseless2DParabola nlp;
    ConjGrad cg(&nlp);

    std::vector<double> initpos{-1., -1.};
    cg.setX(initpos);

    cg.findMin();

    cout << "The found minimum is: ";
    cout << cg.getX(0) << "    " << cg.getX(1) << endl << endl << endl;


    cout << "Now we repeat the minimisation adding a noise to the function and its gradient." << endl;

    cout << "First we try a small amount of noise:" << endl;
    Noisy2DParabola np1(0.05); // sigma = 0.05
    ConjGrad cg2(&np1);
    cg2.setX(initpos);

    cg2.findMin();
    cout << "The found minimum is: ";
    cout << cg2.getX(0) << "    " << cg2.getX(1) << endl << endl;

    cout << "Now we try a larger amount of noise:" << endl;
    Noisy2DParabola np2(0.25); // sigma = 0.25
    ConjGrad cg3(&np2);

    // Let's change some settings for this very noisy CG
    cg3.setStepSize(0.5);
    cg3.setBackStep(0.2); // allow small backstep
    cg3.setX(initpos);

    cg3.findMin();
    cout << "The found minimum is: ";
    cout << cg3.getX(0) << "    " << cg3.getX(1) << endl << endl;

    cout << endl << "As we can see, the CG method still converges to the true minimum,";
    cout << endl << "within the error due to noise. But the method will not improve the";
    cout << endl << "result any further. So in many cases, regular stochastic optimizers";
    cout << endl << "will achieve a better final result (but not in 2 steps!)." << endl;

    // end
    return 0;
}
