#ifndef NFM_NOISYFUNCTION_HPP
#define NFM_NOISYFUNCTION_HPP

#include "nfm/NoisyValue.hpp"
#include "nfm/NoisyGradient.hpp"

#include <vector>

namespace nfm
{

struct NoisyIOPair
// Use this if you want to store pairs
// of input and output of NoisyFunctions.
{
    std::vector<double> x{};
    NoisyValue f{};

    explicit NoisyIOPair(int ndim)
    {
        x.assign(static_cast<size_t>(ndim), 0.);
        f.set(0., 0.);
    }

    int getNDim() const { return static_cast<int>(x.size()); }
};


class NoisyFunction
{
protected:
    const int _ndim;

public:
    explicit NoisyFunction(int ndim): _ndim(ndim) {}
    virtual ~NoisyFunction() = default;

    int getNDim() const { return _ndim; }

    // Noisy Function
    virtual NoisyValue f(const std::vector<double> &x) = 0;
    //         ^ output value&error pair         ^input(size=_ndim)

    // operator () overload
    NoisyValue operator()(const std::vector<double> &x) { return this->f(x); }
};


class NoisyFunctionWithGradient: public NoisyFunction
{
protected:
    const bool _flag_gradErr; // will the function provide gradient errors?

public:
    explicit NoisyFunctionWithGradient(int ndim, bool flag_gradErr):
            NoisyFunction(ndim), _flag_gradErr(flag_gradErr) {}

    bool hasGradErr() const { return _flag_gradErr; }

    // Gradient
    // IMPORTANT: Because we are minimizing, we expect the NEGATIVE gradient.
    virtual void grad(const std::vector<double> &x, NoisyGradient &gradv) = 0;
    //                                         ^input             ^gradient (please set error fields if _flag_gradErr!)

    // Combined Function & Gradient
    // Overwrite it with a more efficient version, if possible
    virtual NoisyValue fgrad(const std::vector<double> &x, NoisyGradient &gradv)
    { //         ^function value&error                ^input             ^gradient output (size ndim)
        NoisyValue ret = this->f(x);
        this->grad(x, gradv);
        return ret;
    }

    NoisyValue operator()(const std::vector<double> &x, NoisyGradient &gradv) { return this->fgrad(x, gradv); }
};
} // namespace nfm

#endif
