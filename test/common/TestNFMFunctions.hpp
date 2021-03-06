#ifndef NFM_TESTNFMFUNCTIONS_HPP
#define NFM_TESTNFMFUNCTIONS_HPP

#include "nfm/NoisyFunction.hpp"
#include <cmath>


class Parabola: public nfm::NoisyFunction
{
public:
    Parabola(): nfm::NoisyFunction(1) {}

    nfm::NoisyValue f(const std::vector<double> &in) override
    {
        return nfm::NoisyValue({in[0]*in[0], 0.5});
    }
};


class Well: public nfm::NoisyFunction
{
public:
    Well(): nfm::NoisyFunction(1) {}

    nfm::NoisyValue f(const std::vector<double> &in) override
    {
        nfm::NoisyValue y{-1, 0.1};
        if ((in[0] <= -1.) || (in[0] >= 1.)) {
            y.val = 1.;
        }
        return y;
    }
};


class PowerFour: public nfm::NoisyFunction
{
public:
    PowerFour(): nfm::NoisyFunction(1) {}

    nfm::NoisyValue f(const std::vector<double> &in) override
    {
        nfm::NoisyValue y{in[0]*in[0], 0.000000001};
        y.val *= y.val;
        return y;
    }
};


class F1D: public nfm::NoisyFunction
{
public:
    F1D(): nfm::NoisyFunction(1) {}

    nfm::NoisyValue f(const std::vector<double> &in) override
    {
        nfm::NoisyValue y{in[0] - 1., 0.00001};
        y.val = pow(y.val, 4);
        return y;
    }
};


class F3D: public nfm::NoisyFunctionWithGradient
{
public:
    F3D(): NoisyFunctionWithGradient(3, true) {}

    nfm::NoisyValue f(const std::vector<double> &in) override   // f = (x-1)^4 + (y+1.5)^4 + (z-0.5)^4
    {
        nfm::NoisyValue y{in[0] - 1., 0.00001};
        y.val = pow(y.val, 4) + pow(in[1] + 1.5, 4) + pow(in[2] - 0.5, 4);
        return y;
    }

    void grad(const std::vector<double> &in, nfm::NoisyGradient &grad) override
    {
        grad.set(0, {-4.*pow(in[0] - 1.0, 3), 0.000001});
        grad.set(1, {-4.*pow(in[1] + 1.5, 3), 0.000001});
        grad.set(2, {-4.*pow(in[2] - 0.5, 3), 0.000001});
    }
};

#endif