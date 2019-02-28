#ifndef NFM_DYNAMICDESCENT_HPP
#define NFM_DYNAMICDESCENT_HPP

#include "nfm/NoisyFunMin.hpp"
#include "nfm/NoisyFunction.hpp"

#include <algorithm>
#include <iostream>


class DynamicDescent: public NFM
{
private:
    double _step_size; // an overall step size scaling factor
    double * _inertia; // inertia per parameter
    double * _old_norm_direction;

    void _writeInertiaInLog();

protected:
    // --- Internal methods
    void findNextX(const double * grad);

public:
    explicit DynamicDescent(NoisyFunctionWithGradient * targetfun, const double stepSize = 1., const bool useGradientError = false, const size_t &max_n_const_values = 20): NFM(targetfun, useGradientError, max_n_const_values), _step_size(stepSize)
    {
        _inertia = new double[_ndim];
        _old_norm_direction = new double[_ndim];
        std::fill(_inertia, _inertia+_ndim, 0.);
        std::fill(_old_norm_direction, _old_norm_direction+_ndim, 0.);

        setGradientTargetFun(targetfun);
    }
    ~DynamicDescent() override{
         delete[] _old_norm_direction;
         delete[] _inertia;
    }

    // --- Minimization
    void findMin() override;
};


#endif
