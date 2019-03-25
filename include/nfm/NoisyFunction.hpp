#ifndef NFM_NOISYFUNCTION_HPP
#define NFM_NOISYFUNCTION_HPP

namespace nfm
{

class NoisyFunction
{
protected:
    int _ndim;

public:
    explicit NoisyFunction(int ndim) { _ndim = ndim; }
    virtual ~NoisyFunction() = default;

    int getNDim() { return _ndim; }

    // Noisy Function
    virtual void f(const double *, double &, double &) = 0;
    //             ^input(size=_ndim)     ^output   ^error on the output
};


class NoisyFunctionWithGradient: public NoisyFunction
{
public:
    explicit NoisyFunctionWithGradient(int ndim): NoisyFunction(ndim) {}
    ~NoisyFunctionWithGradient() override = default;

    // Gradient
    virtual void grad(const double *, double *, double *) = 0;
    //                ^input          ^output   ^error on the output

    // Combined Function & Gradient
    // Overwrite it with a more efficient version, if possible
    virtual void fgrad(const double * input, double &output, double &error, double * gradient, double * graderr)
    {
        f(input, output, error);
        grad(input, gradient, graderr);
    }
};


//class NoisyFunctionWithDomain: public NoisyFunction
//{
//   public:
//      // Domain
//      virtual bool domain(const double *) = 0;
//      //                  ^input(size=_ndim)
//};
} // namespace nfm

#endif
