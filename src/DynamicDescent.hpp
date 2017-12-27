#ifndef DYNAMIC_DESCENT
#define DYNAMIC_DESCENT

#include "NoisyFunMin.hpp"
#include "NoisyFunction.hpp"

class DynamicDescent: public NFM{
   
private:
   double _inertia;
   double * _old_norm_direction;
   
   void writeCurrentXInLog();
   void writeDirectionInLog(const double * direction, const double * directionerror);
   void reportMeaninglessGradientInLog();
   void writeInertiaInLog();
   
protected:

   // --- Internal methods
   void findNextX(const double * dir, double &deltatargetfun, double &deltax);

public:
   DynamicDescent(NoisyFunctionWithGradient * targetfun):NFM(targetfun)
   {
      _inertia = 1.;
      _old_norm_direction = new double[targetfun->getNDim()];
      for (int i=0; i<targetfun->getNDim(); ++i){
         _old_norm_direction[i] = 0.;
      }
      setGradientTargetFun(targetfun);
   }
   ~DynamicDescent(){
      delete[] _old_norm_direction;
   }

   // --- Minimization
   void findMin();

};


#endif
