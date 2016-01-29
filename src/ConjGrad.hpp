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
      ConjGrad(NoisyFunctionWithGradient * targetfun):NFM(targetfun)
      {
         setGradientTargetFun(targetfun);
      }
      ~ConjGrad(){ }


      // --- Minimization
      void findMin();
      
};


#endif
