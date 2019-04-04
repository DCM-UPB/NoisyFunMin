#ifndef NFM_DYNAMICDESCENT_HPP
#define NFM_DYNAMICDESCENT_HPP

#include "nfm/NoisyFunMin.hpp"

namespace nfm
{

enum class DDMode
{
    SGDM, /* SGD with momentum */
    ADAG, /* AdaGrad */
    ADAD, /* AdaDelta */
    RMSP, /* RMSProp */
    NEST  /* Nesterov */
};

// Stochastic Gradient Descent Algorithms
//
// All contained algorithms provide some kind of adaptive learning rate,
// usually by the use of momenta. There
class DynamicDescent: public NFM
{
private:
    DDMode _ddmode; // which update rule to use
    bool _useAveraging; // use the averaged positions of the old value list (length max_n_const_values) as end result
    double _stepSize; // step size factor / learning rate
    double _beta; // momenta update parameter (not used in AdaGrad)
    double _epsilon; // small value to prevent bad division (not used in SGDM)

    // --- Internal methods
    bool _updateTarget();
    void _findNextX(int iter, std::vector<double> &v, std::vector<double> &w);
    void _findMin() override;

public:
    explicit DynamicDescent(NoisyFunctionWithGradient * targetfun, DDMode ddmode = DDMode::SGDM, bool useAveraging = false, double stepSize = 0.01, double beta = 0.9, double epsilon = 1.e-8);
    ~DynamicDescent() override = default;

    // DD Configuration
    void useSGDM() { _ddmode = DDMode::SGDM; }
    void useAdaGrad() { _ddmode = DDMode::ADAG; }
    void useAdaDelta() { _ddmode = DDMode::ADAD; }
    void useRMSProp() { _ddmode = DDMode::RMSP; }
    void useNesterov() { _ddmode = DDMode::NEST; }
    void setDDMode(DDMode ddmode) { _ddmode = ddmode; }
    DDMode getDDMode() const { return _ddmode; }

    // Getters
    bool usesAveraging() const { return _useAveraging; }
    double getStepSize() const { return _stepSize; }
    double getBeta() const { return _beta; }

    // Setters
    void setAveraging(bool useAveraging) { _useAveraging = useAveraging; }
    void setStepSize(double stepSize) { _stepSize = std::max(0., stepSize); }
    void setBeta(double beta) { _beta = std::max(0., std::min(1., beta)); }
};
} // namespace nfm

#endif
