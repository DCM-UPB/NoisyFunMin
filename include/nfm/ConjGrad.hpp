#ifndef NFM_CONJGRAD_HPP
#define NFM_CONJGRAD_HPP

#include "nfm/NoisyFunMin.hpp"
#include "nfm/NoisyFunction.hpp"

namespace nfm
{

enum class CGMode
{
    RAW /* use raw gradients */,
    CGFR /* use Fletcher-Reeves CG */
};

// Noisy Conjugate-Gradient Minimization
class ConjGrad: public NFM
{
private:
    double _stepSize; // search interval size
    CGMode _cgmode; // which gradients to use

    // --- Internal methods
    void _findNextX(const std::vector<double> &dir); // do line-search
    void _writeCGDirectionToLog(const std::vector<double> &dir, const std::string &name) const;

    // --- Minimization
    void _findMin() final; // perform noisy CG minimization

public:
    explicit ConjGrad(NoisyFunctionWithGradient * targetfun, int max_n_const_values = 20, double stepSize = 0.):
            NFM(targetfun, max_n_const_values), _stepSize(stepSize), _cgmode(CGMode::CGFR) {}
    ~ConjGrad() final = default;


    // Configuration
    void useRawGrad() { _cgmode = CGMode::RAW; }  // make ConjGrad Steepest-Descent-like
    void useConjGradFR() { _cgmode = CGMode::CGFR; }  // reset to normal
    void setCGMode(CGMode cgmode) { _cgmode = cgmode; }
    CGMode getCGMode() const { return _cgmode; }
};
} // namespace nfm

#endif
