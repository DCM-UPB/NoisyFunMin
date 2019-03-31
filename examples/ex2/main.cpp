#include "nfm/DynamicDescent.hpp"
#include "nfm/LogManager.hpp"

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
    //LogManager::setLoggingOn(true); // use this for verbose printout of the DD method

    cout << "We first minimize it, supposing to have no noise at all" << endl;

    Noiseless3DParabola nlp;
    DynamicDescent dd(&nlp);

    double initpos[3]{2.5, 1., -1.};
    dd.setX(initpos);

    // in the case without noise we change settings
    dd.setStepSize(0.5); // a larger step size
    dd.setEpsF(0.001); // we should enable stopping on too small function changes

    // and now find the min
    dd.findMin();

    cout << "The found minimum is: ";
    cout << dd.getX(0) << "    " << dd.getX(1) << "    " << dd.getX(2) << endl << endl;


    cout << "Now we repeat the minimisation adding a noise to the function and its gradient." << endl;

    NoisyWrapper np(&nlp, 0.25); // sigma 0.25
    DynamicDescent dd2(&np);
    dd2.setStepSize(0.01); // a smaller step size
    dd2.setMaxNConstValues(20);
    dd2.setAveraging(true); // this option usually improves the final result significantly
    dd2.setX(initpos);

    dd2.findMin();

    cout << "The found minimum is: ";
    cout << dd2.getX(0) << "    " << dd2.getX(1) << "    " << dd2.getX(2) << endl << endl;


    // end
    return 0;
}
