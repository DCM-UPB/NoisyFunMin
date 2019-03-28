#include <cassert>
#include <cmath>
#include <iostream>

#include "nfm/Adam.hpp"
#include "nfm/LogManager.hpp"

#include "TestNFMFunctions.hpp"


int main()
{
    using namespace std;
    using namespace nfm;

    //LogManager::setLogLevel(LogLevel::VERBOSE);

    // define 3D function that I want to minimise
    F3D f3d;

    // introduce array with the initial position
    double x[3];

    bool useGradientError = false;
    bool useAveraging = true;
    for (int i = 0; i < 2; ++i) {
        useGradientError = !useGradientError;
        for (int j = 0; j < 2; ++j) {
            useAveraging = !useAveraging;

            // cout << "useGradientError " << useGradientError << " useAveraging " << useAveraging << endl;

            // test Adam
            Adam adam = Adam(&f3d, 20, useAveraging, 0.1, 0.1, 0.1); // since the gradient is exact we chose very high decay (0 would actually be ideal)
            x[0] = -2.;
            x[1] = 1.0;
            x[2] = 0.0;
            adam.setX(x);
            adam.findMin();

            assert(fabs(adam.getX(0) - 1.0) < 0.1);
            assert(fabs(adam.getX(1) + 1.5) < 0.1);
            assert(fabs(adam.getX(2) - 0.5) < 0.1);
        }
    }


    return 0;
}
