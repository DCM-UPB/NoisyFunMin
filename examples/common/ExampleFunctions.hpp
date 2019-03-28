#ifndef NFM_EXAMPLEFUNCTIONS_HPP
#define NFM_EXAMPLEFUNCTIONS_HPP

#include "nfm/NoisyFunction.hpp"

#include <cmath>
#include <random>


class Noiseless2DParabola: public nfm::NoisyFunctionWithGradient
{
public:
    Noiseless2DParabola(): nfm::NoisyFunctionWithGradient(2, true) {}

    nfm::NoisyValue f(const std::vector<double> &in) override
    {
        nfm::NoisyValue y{};
        y.value = pow(in[0] - 1., 2) + pow(in[1] + 2., 2);   // minimum in (1, -2)
        return y;
    }

    void grad(const std::vector<double> &in, std::vector<nfm::NoisyValue> &grad) override
    {
        grad[0].value = 2.*(in[0] - 1.);
        grad[1].value = 2.*(in[1] + 2.);
    }
};


class Noisy2DParabola: public nfm::NoisyFunctionWithGradient
{
private:
    const double _sigma = 0.25;
    std::random_device _rdev;
    std::mt19937_64 _rgen;
    std::uniform_real_distribution<double> _rd;  //after initialization (done in the constructor) can be used with _rd(_rgen)

public:
    Noisy2DParabola(): nfm::NoisyFunctionWithGradient(2, true)
    {
        // initialize random generator
        _rgen = std::mt19937_64(_rdev());
        _rd = std::uniform_real_distribution<double>(-_sigma, _sigma);
    }

    nfm::NoisyValue f(const std::vector<double> &in) override
    {
        nfm::NoisyValue y{_rd(_rgen), _sigma};
        y.value += pow(in[0] - 1., 2) + pow(in[1] + 2., 2);   // minimum in (1, -2)
        return y;
    }

    void grad(const std::vector<double> &in, std::vector<nfm::NoisyValue> &grad) override
    {
        grad[0].value = 2.*(in[0] - 1.) + 2.*_rd(_rgen);
        grad[1].value = 2.*(in[1] + 2.) + 2.*_rd(_rgen);
        grad[0].error = 2.*_sigma;
        grad[1].error = 2.*_sigma;
    }
};


#endif