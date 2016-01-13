#ifndef CONJ_GRAD
#define CONJ_GRAD

#include "NoisyFunMin.hpp"
#include "NoisyFunction.hpp"

class ConjGrad: public NFM
{
   protected:
      
      // --- Internal methods
      void findNextX(const double * dir, double &deltatargetfun, double &deltax);

   public:
      ConjGrad(const int &ndim, NoisyFunctionWithGradient * targetfun):NFM(ndim,targetfun)
      {
         setGradientTargetFun(targetfun);
      }
      ~ConjGrad(){ }


      // --- Minimization
      void findMin();
      
};


#endif
