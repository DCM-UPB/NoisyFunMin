#ifndef NOISY_FUNCTION
#define NOISY_FUNCTION


class NoisyFunction
{
protected:
    int _ndim;

public:
    NoisyFunction(int ndim){_ndim=ndim;}
    virtual ~NoisyFunction(){}

    int getNDim(){return _ndim;}

    // Noisy Function
    virtual void f(const double *, double &, double &) = 0;
    //             ^input(size=_ndim)     ^output   ^error on the output
};


class NoisyFunctionWithGradient: public NoisyFunction
{
public:
    NoisyFunctionWithGradient(int ndim): NoisyFunction(ndim){}
    virtual ~NoisyFunctionWithGradient(){}

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


#endif
