#include "nfm/Adam.hpp"
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
    cout << "Adam Example" << endl << endl;
    cout << "We want to minimize the 3D function" << endl;
    cout << "    x^2 + (y+1)^2 + (z-2)^2" << endl;
    cout << "whose min is in (0, -1, 2)." << endl << endl;
    cout << "We will always start at (2.5, 1, -1)." << endl << endl << endl;

    cout << "We first minimize it, supposing to have no noise at all" << endl;

    Noiseless3DParabola nlp;
    Adam adam(&nlp);

    // Make sure that the noiseless case converges quickly
    // Note that using Adam for noiseless optimization is abuse
    adam.setAlpha(0.5);
    adam.setBeta1(0.1);
    adam.setBeta2(0.5);
    adam.setEpsF(0.01); // stop on too small function changes

    adam.setX({2.5, 1., -1.});
    adam.findMin();
    reportMinimum(adam);


    cout << "Now we repeat the minimisation adding a noise to the function and its gradient." << endl;

    NoisyWrapper np(&nlp, 0.25); // sigma 0.25
    Adam adam2(&np, true /* use averaging to calculate final parameters */, 0.5 /* step size factor */);

    adam2.setX(adam.getX());
    adam2.findMin();
    reportMinimum(adam2);

    cout << "NOTE: You may enable detailed logging by uncommenting" << endl;
    cout << "      one of two lines in the beginning of the example." << endl;

    // end
    return 0;
}

