#ifndef DYNAMIC_DESCENT
#define DYNAMIC_DESCENT

#include "nfm/NoisyFunMin.hpp"
#include "nfm/NoisyFunction.hpp"
#include "nfm/NoisyFunctionValue.hpp"

#include <list>
#include <iostream>


class DynamicDescent: public NFM{

private:
    double _step_size; // an overall step size scaling factor
    double * _inertia; // inertia per parameter
    double * _old_norm_direction;

    void _writeInertiaInLog();

protected:
    // --- Internal methods
    void findNextX(const double * dir);

public:
    DynamicDescent(NoisyFunctionWithGradient * targetfun, const double stepSize = 1., const bool useGradientError = false, const size_t &max_n_const_values = 20): NFM(targetfun, useGradientError, max_n_const_values), _step_size(stepSize)
    {
        _inertia = new double[_ndim];
        _old_norm_direction = new double[_ndim];
        for (int i=0; i<_ndim; ++i){
            _inertia[i] = 0.;
            _old_norm_direction[i] = 0.;
        }
        setGradientTargetFun(targetfun);
    }
    ~DynamicDescent(){ delete[] _old_norm_direction; }

    // --- Minimization
    void findMin();
};


#endif
