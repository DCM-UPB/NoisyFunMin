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
    LogManager::setLogLevel(LogLevel::VERBOSE);

    // define 3D function that I want to minimise
    F3D f3d;

    for (int i = 0; i < 2; ++i) {
        const bool useGradientError = (i != 0); // this option isn't Adam-specific
        for (int j = 0; j < 2; ++j) {
            const bool useAveraging = (j != 0);
            for (int k = 0; k < 2; ++k) {
                const bool useAMSGrad = (k != 0);

                //cout << "useGradientError " << useGradientError << " useAveraging " << useAveraging << " useAMSGrad " << useAMSGrad << endl;

                // test Adam
                Adam adam(f3d.getNDim(), useAveraging, 0.1);
                adam.setMaxNConstValues(100);
                if (!useAMSGrad) {
                    adam.setBeta1(0.1); // this case seems to work better with high decay
                    adam.setBeta2(0.1);
                }
                else {
                    adam.setAlpha(1.);
                    adam.setBeta1(0.5);
                    adam.setBeta2(0.9);
                }
                adam.setGradErrStop(useGradientError);
                adam.setAMSGrad(useAMSGrad);
                adam.setX(0, -2.);
                adam.setX(1, 1.);
                adam.setX(2, 0.);
                NoisyIOPair opt = adam.findMin(f3d);

                assert(fabs(opt.x[0] - 1.0) < 0.1);
                assert(fabs(opt.x[1] + 1.5) < 0.1);
                assert(fabs(opt.x[2] - 0.5) < 0.1);
            }
        }
    }

    return 0;
}
