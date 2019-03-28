#include "nfm/DynamicDescent.hpp"
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

    //LogManager::setLoggingOn(); // use this to enable log printout
    //LogManager::setLoggingOn(true); // use this for verbose printout of the DD method

    cout << "we first minimize it, supposing to have no noise at all" << endl;

    Noiseless2DParabola nlp;
    DynamicDescent dd(&nlp);

    double initpos[2] {-1., 1.};
    dd.setX(initpos);

    dd.findMin();

    cout << "The found minimum is: ";
    cout << dd.getX(0) << "    " << dd.getX(1) << endl << endl << endl;


    cout << "Now we repeat the minimisation adding a noise to the function and its gradient." << endl;

    Noisy2DParabola np;
    DynamicDescent dd2(&np);

    dd2.setX(initpos);

    dd2.findMin();

    cout << "The found minimum is: ";
    cout << dd2.getX(0) << "    " << dd2.getX(1) << endl << endl;


    // end
    return 0;
}
