#ifndef DYNAMIC_DESCENT
#define DYNAMIC_DESCENT

#include "NoisyFunMin.hpp"
#include "NoisyFunction.hpp"
#include "NoisyFunctionValue.hpp"

#include <list>
#include <iostream>


class DynamicDescent: public NFM{

private:
    double _inertia;
    double * _old_norm_direction;
    std::list<NoisyFunctionValue *> _old_values;
    const unsigned int N_CONSTANT_VALUES_CONDITION_FOR_STOP = 20;

    void writeCurrentXInLog();
    void writeDirectionInLog(const double * direction, const double * directionerror);
    void reportMeaninglessGradientInLog();
    void writeInertiaInLog();
    void writeOldValuesInLog();

protected:

    // --- Internal methods
    void findNextX(const double * dir);
    bool shouldContinueDescent();

public:
    DynamicDescent(NoisyFunctionWithGradient * targetfun):NFM(targetfun)
    {
        _inertia = 1.;
        _old_norm_direction = new double[targetfun->getNDim()];
        for (int i=0; i<targetfun->getNDim(); ++i){
            _old_norm_direction[i] = 0.;
        }
        _old_values.clear();
        setGradientTargetFun(targetfun);
    }
    ~DynamicDescent(){
        for (NoisyFunctionValue * v : _old_values)
            delete v;
        _old_values.clear();
        delete[] _old_norm_direction;
    }

    // --- Minimization
    void findMin();

};


#endif
