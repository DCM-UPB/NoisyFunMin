#ifndef NFM_EXAMPLEFUNCTIONS_HPP
#define NFM_EXAMPLEFUNCTIONS_HPP

#include "nfm/NoisyFunction.hpp"

#include <cmath>
#include <random>


class Noiseless3DParabola: public nfm::NoisyFunctionWithGradient
{
public:
    Noiseless3DParabola(): nfm::NoisyFunctionWithGradient(3, true) {}

    nfm::NoisyValue f(const std::vector<double> &in) override
    {
        nfm::NoisyValue y{};
        y.value = pow(in[0], 2) + pow(in[1] + 1., 2) + pow(in[2] - 2., 2);   // minimum in (0, -1, 2)
        return y;
    }

    void grad(const std::vector<double> &in, std::vector<nfm::NoisyValue> &grad) override
    {
        grad[0].value = 2*in[0];
        grad[1].value = 2.*(in[1] + 1.);
        grad[2].value = 2.*(in[2] - 2.);
    }
};

class NoisyWrapper: public nfm::NoisyFunctionWithGradient
{
private:
    nfm::NoisyFunctionWithGradient * const _nlfun; // the noiseless function to wrap
    double _sigma;
    std::random_device _rdev;
    std::mt19937_64 _rgen;
    std::normal_distribution<double> _rd; // this returns random doubles when called as_rd(_rgen)

public:
    explicit NoisyWrapper(nfm::NoisyFunctionWithGradient * fun, double sigma = 1.):
    nfm::NoisyFunctionWithGradient(fun->getNDim(), fun->hasGradErr()), _nlfun(fun), _sigma(sigma)
    {
        // initialize random generator
        _rgen = std::mt19937_64(_rdev());
        _rd = std::normal_distribution<double>(0., _sigma);
    }

    nfm::NoisyValue makeValueNoisy(nfm::NoisyValue nv, double sigfac = 1.) // a factor to increase error
    {
        nv.value += sigfac*_rd(_rgen);
        nv.error = sigfac*_sigma;
        return nv;
    }

    nfm::NoisyValue f(const std::vector<double> &in) final
    {
        return makeValueNoisy(_nlfun->f(in)); // use helper function above
    }

    void grad(const std::vector<double> &in, std::vector<nfm::NoisyValue> &grad) final
    {
        _nlfun->grad(in, grad);
        for (auto & gi : grad) { gi = makeValueNoisy(gi, 2.); } // gradients have larger error
    }
};


#endif
