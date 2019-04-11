#include <cassert>
#include <cmath>

#include "nfm/Adam.hpp"
#include "nfm/LogManager.hpp"

#include "TestNFMFunctions.hpp"


int main()
{
    using namespace std;
    using namespace nfm;

    LogManager::setLoggingOff();
    //LogManager::setLogLevel(LogLevel::VERBOSE);

    // define 3D function that I want to minimise
    F3D f3d;

    for (int i = 0; i < 2; ++i) {
        const bool useGradientError = (i != 0); // this option isn't Adam-specific
        for (int j = 0; j < 2; ++j) {
            const bool useAveraging = (j != 0);

            // cout << "useGradientError " << useGradientError << " useAveraging " << useAveraging << endl;

            // test Adam
            Adam adam(&f3d, useAveraging, 0.1);
            adam.setBeta1(0.1); // this case works better with high decay
            adam.setBeta2(0.1);
            adam.setGradErrStop(useGradientError);
            adam.setX(0, -2.);
            adam.setX(1, 1.);
            adam.setX(2, 0.);
            adam.findMin();

            assert(fabs(adam.getX(0) - 1.0) < 0.1);
            assert(fabs(adam.getX(1) + 1.5) < 0.1);
            assert(fabs(adam.getX(2) - 0.5) < 0.1);
        }
    }

    return 0;
}
