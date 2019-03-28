#ifndef NFM_DYNAMICDESCENT_HPP
#define NFM_DYNAMICDESCENT_HPP

#include "nfm/NoisyFunMin.hpp"
#include "nfm/NoisyFunction.hpp"

#include <algorithm>
#include <iostream>
#include <cmath>

namespace nfm
{

class DynamicDescent: public NFM
{
private:
    double _stepSize; // an overall step size scaling factor
    std::vector<double> _inertia; // inertia per parameter
    std::vector<double> _old_norm_dir;

    void _writeInertiaToLog();

protected:
    // --- Internal methods
    void findNextX(const std::vector<NoisyValue> &grad);
    void _findMin() override;

public:
    explicit DynamicDescent(NoisyFunctionWithGradient * targetfun, double stepSize = 0.01, int max_n_const_values = 20);

    double getStepSize() const { return _stepSize; }
    void setStepSize(double stepSize) { _stepSize = stepSize; }
};
} // namespace nfm

#endif
