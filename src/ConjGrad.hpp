#ifndef CONJ_GRAD
#define CONJ_GRAD

#include "NoisyFunMin.hpp"
#include "NoisyFunction.hpp"

class ConjGrad: public NFM{

private:
    bool _use_conjgrad;

    void writeCurrentXInLog();
    void writeDirectionInLog(const double * direction, const double * directionerror);
    void reportMeaninglessGradientInLog();

protected:

    // --- Internal methods
    void findNextX(const double * dir, double &deltatargetfun, double &deltax);

public:
    ConjGrad(NoisyFunctionWithGradient * targetfun):NFM(targetfun)
    {
        setGradientTargetFun(targetfun);
        _use_conjgrad = true;
    }
    ~ConjGrad(){ }


    // Configuration
    void configureToFollowSimpleGradient(){_use_conjgrad = false;};  // make ConjGrad a Steepest Descent


    // --- Minimization
    void findMin();

};


#endif
