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

    Noisy2DParabola np;
    ConjGrad cg2(&np);
    cg2.setX(initpos);

    // For noisy CG we might want to increase
    // the tolerances / decrease target precision.
    cg2.setEpsX(0.01);
    cg2.setEpsF(0.01);

    cg2.findMin();

    cout << "The found minimum is: ";
    cout << cg2.getX(0) << "    " << cg2.getX(1) << endl << endl;

    // end
    return 0;
}
