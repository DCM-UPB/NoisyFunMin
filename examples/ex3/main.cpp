#include "nfm/Adam.hpp"
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
    //LogManager::setLoggingOn(true); // use this for verbose printout of the Adam method

    cout << "We first minimize it, supposing to have no noise at all" << endl;

    Noiseless3DParabola nlp;

    Adam adam(&nlp);
    adam.setX({2.5, 1., -1.});

    // Make sure that the noiseless case converges quickly
    // Note that using Adam for noiseless optimization is abuse
    adam.setAlpha(0.5);
    adam.setBeta1(0.1);
    adam.setBeta2(0.5);
    adam.setEpsF(0.01); // stop on too small function changes
    adam.findMin();

    cout << "The found minimum is: ";
    cout << adam.getX(0) << "    " << adam.getX(1) << "    " << adam.getX(2) << endl << endl;


    cout << "Now we repeat the minimisation adding a noise to the function and its gradient." << endl;

    NoisyWrapper np(&nlp, 0.25); // sigma 0.25

    Adam adam2(&np, true /* use averaging to calculate final parameters */, 0.5 /* step size factor */);
    adam2.setX({2.5, 1., -1.});
    adam2.findMin();

    cout << "The found minimum is: ";
    cout << adam2.getX(0) << "    " << adam2.getX(1) << "    " << adam2.getX(2) << endl << endl;


    // end
    return 0;
}

