#ifndef ADAM
#define ADAM

#include "NoisyFunMin.hpp"
#include "NoisyFunction.hpp"
#include "NoisyFunctionValue.hpp"

#include <list>

class Adam: public NFM
{
protected:
    double _alpha; // stepsize, default 0.001
    double _beta1, _beta2; // decay rates in [0, 1), default 0.9 and 0.999 respectively
    double _epsilon; // offset to stabilize division in update, default 10e-8

public:
    Adam(NoisyFunctionWithGradient * targetfun, const bool useGradientError = false, const size_t &max_n_const_values = 20, const double &alpha = 0.001, const double &beta1 = 0.9, const double &beta2 = 0.999, const double &epsilon = 10e-8)
        : NFM(targetfun, useGradientError, max_n_const_values), _alpha(alpha), _beta1(beta1), _beta2(beta2), _epsilon(epsilon) { setGradientTargetFun(targetfun); }

    // --- Minimization
    void findMin();
};


#endif
