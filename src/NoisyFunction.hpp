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
};


//class NoisyFunctionWithDomain: public NoisyFunction
//{
//   public:
//      // Domain
//      virtual bool domain(const double *) = 0;
//      //                  ^input(size=_ndim)
//};


#endif
