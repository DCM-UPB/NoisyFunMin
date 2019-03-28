#include "nfm/Adam.hpp"
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
    //LogManager::setLoggingOn(true); // use this for verbose printout of the Adam method

    cout << "we first minimize it, supposing to have no noise at all" << endl;

    Noiseless2DParabola nlp;

    Adam adam(&nlp);

    adam.setX(0, -1.);
    adam.setX(1, -1.);

    // make sure that the noiseless case converges quickly
    adam.setAlpha(0.1);
    adam.setBeta1(0.1);
    adam.setBeta2(0.5);
    adam.setEpsF(0.001);
    adam.findMin();

    cout << "The found minimum is: ";
    cout << adam.getX(0) << "    " << adam.getX(1) << endl << endl << endl;


    cout << "Now we repeat the minimisation adding a noise to the function and its gradient." << endl;

    Noisy2DParabola np;

    // in noisy-target low-dim cases like this we often need to change adam parameters from default (in brackets)
    Adam adam2(&np, 20 /* max n constant (within error) values before stopping */,
               true /* use averaging to calculate final parameters */, 1.0 /* step size factor (0.001) */);
    adam2.setX({-1., -1.});

    adam2.findMin();

    cout << "The found minimum is: ";
    cout << adam2.getX(0) << "    " << adam2.getX(1) << endl << endl;


    // end
    return 0;
}
