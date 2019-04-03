#ifndef NFM_EXAMPLEFUNCTIONS_HPP
#define NFM_EXAMPLEFUNCTIONS_HPP

#include "nfm/NoisyFunction.hpp"

#include <cmath>
#include <random>

// --- The Noisy Wrapper

class NoisyWrapper: public nfm::NoisyFunctionWithGradient
{
private:
    nfm::NoisyFunctionWithGradient * const _nlfun; // the noiseless function to wrap
    std::random_device _rdev;
    std::mt19937_64 _rgen;
    std::normal_distribution<double> _rd; // this returns random doubles when called as_rd(_rgen)
    double _sigma; // the standard deviation

public:
    explicit NoisyWrapper(nfm::NoisyFunctionWithGradient * fun, double sigma_noise = 1.):
            nfm::NoisyFunctionWithGradient(fun->getNDim(), true), _nlfun(fun), _sigma(sigma_noise)
    {
        // initialize random generator
        _rgen = std::mt19937_64(_rdev());
        _rd = std::normal_distribution<double>(0., _sigma);
    }

    double getSigma() const { return _sigma; }
    void setSigma(double sigma)
    {
        _sigma = sigma;
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
        for (auto &gi : grad) { gi = makeValueNoisy(gi, 2.); } // gradients have larger error
    }
};


// --- Noiseless Test Functions

class TestParabola3D: public nfm::NoisyFunctionWithGradient
{
public:
    TestParabola3D(): nfm::NoisyFunctionWithGradient(3, false/*no grad errors*/) {}

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


template<int N>
class RosenbrockFunction: public nfm::NoisyFunctionWithGradient
{
public:
    RosenbrockFunction(): nfm::NoisyFunctionWithGradient(N, false) {}

    nfm::NoisyValue f(const std::vector<double> &in) override
    {
        double y = 0;
        for (int i=0; i<N-1; ++i) { y += 100.*pow(in[i+1]-in[i]*in[i], 2) + pow(1.-in[i], 2);}
        return {y, 0.};
    }

    void grad(const std::vector<double> &in, std::vector<nfm::NoisyValue> &grad) override
    {
        std::fill(grad.begin(), grad.end(), nfm::NoisyValue({0., 0.}));
        for (int i=0; i<N-1; ++i) {
            const double common = 200.*(in[i+1]-in[i]*in[i]);
            grad[i].value -= common*2.*in[i] + 2.*(1.-in[i]);
            grad[i+1].value += common;
        }
    }
};

#endif
