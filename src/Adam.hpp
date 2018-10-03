#ifndef ADAM
#define ADAM

#include "NoisyFunMin.hpp"
#include "NoisyFunction.hpp"
#include "NoisyFunctionValue.hpp"

#include <list>

class Adam: public NFM
{
protected:
    double _alpha; // stepsize, default 0.001
    double _beta1, _beta2; // decay rates in [0, 1), default 0.9 and 0.999 respectively
    double _epsilon; // offset to stabilize division in update, default 10e-8

    std::list<NoisyFunctionValue *> _old_values;
    const unsigned int N_CONSTANT_VALUES_CONDITION_FOR_STOP = 20;

    void writeCurrentXInLog();
    void writeDirectionInLog(const double * direction);
    void reportMeaninglessGradientInLog();

    // --- Internal methods
    bool isNotConverged();

public:
    Adam(NoisyFunctionWithGradient * targetfun, const double &alpha = 0.001, const double &beta1 = 0.9, const double &beta2 = 0.999, const double &epsilon = 10e-8)
        : NFM(targetfun), _alpha(alpha), _beta1(beta1), _beta2(beta2), _epsilon(epsilon) { setGradientTargetFun(targetfun); }

    ~Adam()
    {
        for (NoisyFunctionValue * v : _old_values) {
            delete v;
        }
        _old_values.clear();
    }

    // --- Minimization
    void findMin();

};


#endif
