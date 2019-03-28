#ifndef NFM_DYNAMICDESCENT_HPP
#define NFM_DYNAMICDESCENT_HPP

#include "nfm/NoisyFunMin.hpp"
#include "nfm/NoisyFunction.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>

namespace nfm
{

// Basic Stochastic Gradient Descent with Momenta
class DynamicDescent: public NFM
{
private:
    double _stepSize; // step size factor / learning rate
    double _alpha; // momenta update step-size

    // --- Internal methods
    void _findNextX(const std::vector<NoisyValue> &grad, std::vector<double> &dx);
    void _findMin() final;

public:
    explicit DynamicDescent(NoisyFunctionWithGradient * targetfun, int max_n_const_values = DEFAULT_MAX_N_CONST,
                            double stepSize = 0.01, double alpha = 0.9);
    ~DynamicDescent() final = default;

    double getStepSize() const { return _stepSize; }
    void setStepSize(double stepSize) { _stepSize = std::max(0., stepSize); }
    double getAlpha() const { return _alpha; }
    void setAlpha(double alpha) { _alpha = std::max(0., alpha); }
};
} // namespace nfm

#endif
