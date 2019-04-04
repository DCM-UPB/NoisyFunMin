#ifndef NFM_ADAM_HPP
#define NFM_ADAM_HPP

#include "nfm/NoisyFunMin.hpp"

namespace nfm
{

// Adam algorithm, based on https://arxiv.org/abs/1412.6980
//
// This method is similar in principle to the optimization
// methods provided by the DynamicDescent class. However,
// Adam has a more complex update scheme involving second
// order momentum and it provides an inherent averaging method.
class Adam: public NFM
{
private:
    bool _useAveraging; // use automatic exponential decaying (beta2) parameter averaging, as proposed in the end of Adam paper
    double _alpha; // stepsize, default 0.001
    double _beta1 = 0.9, _beta2 = 0.999; // decay rates in [0, 1)
    double _epsilon = 1.e-8; // offset to stabilize division in update

    // --- Minimization
    void _findMin() override;

public:
    explicit Adam(NoisyFunctionWithGradient * targetfun, bool useAveraging = false, double alpha = 0.001);
    ~Adam() override = default;

    // Getters
    bool usesAveraging() const { return _useAveraging; }
    double getAlpha() const { return _alpha; }
    double getBeta1() const { return _beta1; }
    double getBeta2() const { return _beta2; }
    double get_epsilon() const { return _epsilon; }

    // Setters
    void setAveraging(bool useAveraging) { _useAveraging = useAveraging; }
    void setAlpha(double alpha) { _alpha = std::max(0., alpha); }
    void setBeta1(double beta1) { _beta1 = std::max(0., std::min(1., beta1)); }
    void setBeta2(double beta2) { _beta2 = std::max(0., std::min(1., beta2)); }
    void setEpsilon(double epsilon) { _epsilon = std::max(0., epsilon); }
};
} // namespace nfm

#endif
