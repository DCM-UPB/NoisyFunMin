#include "nfm/DynamicDescent.hpp"
#include "nfm/Adam.hpp"
#include "nfm/LogManager.hpp"

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
    cout << "Stochastic Gradient Descent Example" << endl << endl;
    cout << "We want to minimize the 3D function" << endl;
    cout << "    x^2 + (y+1)^2 + (z-2)^2" << endl;
    cout << "whose min is in (0, -1, 2)." << endl << endl;
    cout << "We will always start at (2.5, 1, -1)." << endl << endl << endl;

    cout << "We first minimize it, supposing to have no noise at all" << endl;

    TestParabola3D nlp;
    DynamicDescent dd(&nlp);
    double initpos[3]{2.5, 1., -1.};

    // in the case without noise we change settings
    dd.setStepSize(0.5); // the magic step size
    dd.setBeta(0.); // disable momentum
    dd.setEpsF(0.001); // we should enable stopping on too small function changes

    // and now find the min
    dd.setX(initpos);
    dd.findMin();
    reportMinimum(dd);


    cout << "Now we repeat the minimisation adding a noise to the function and its gradient." << endl;

    NoisyWrapper np(&nlp, 0.25); // sigma 0.25
    DynamicDescent dd2(&np);

    // here we use different settings
    dd2.setStepSize(0.01); // a smaller step size
    dd2.setMaxNConstValues(20);
    dd2.setAveraging(true); // this option usually improves the final result significantly

    dd2.setX(initpos);
    dd2.findMin();
    reportMinimum(dd2);

    cout << "We may also use a different Stochastic Gradient algorithm, like AdaDelta:" << endl;

    dd2.useAdaDelta();
    dd2.setStepSize(0.1); // the initial step size (in AdaDelta this is only used in step 1!)

    dd2.setX(initpos);
    dd2.findMin();
    reportMinimum(dd2);


    cout << "Another SGD algorithm is Adam, which resides in its own class. Let's try it:" << endl;

    Adam adam(&np, true /* use averaging to obtian final result */, 0.1 /* step size factor */);

    adam.setX(dd2.getX());
    adam.findMin();
    reportMinimum(adam);

    cout << "NOTE: You may enable detailed logging by uncommenting" << endl;
    cout << "      one of two lines in the beginning of the example." << endl;

    // end
    return 0;
}
