#ifndef NFM_CONJGRAD_HPP
#define NFM_CONJGRAD_HPP

#include "nfm/NoisyFunMin.hpp"
#include "nfm/NoisyFunction.hpp"

namespace nfm
{

enum class CGMode
{
    NOCG, /* use raw gradients */
    CGFR, /* use Fletcher-Reeves CG */
    CGPR, /* use Polak-Ribiere CG */
    CGPR0  /* use PR-CG with reset */
};

// Noisy Conjugate-Gradient Minimization
class ConjGrad: public NFM
{
private:
    double _stepSize; // initial search interval size factor
    int _max_n_bracketing; // maximal amount of iterations in findBracket()
    CGMode _cgmode; // which gradients to use

    // --- Internal methods
    void _findNextX(const std::vector<double> &dir); // do line-search
    void _writeCGDirectionToLog(const std::vector<double> &dir, const std::string &name) const;

    // --- Minimization
    void _findMin() final; // perform noisy CG minimization

public:
    explicit ConjGrad(NoisyFunctionWithGradient * targetfun, double stepSize = 1., int max_n_bracketing = 10);
    ~ConjGrad() final = default;

    // CG Configuration
    void useRawGrad() { _cgmode = CGMode::NOCG; }  // make ConjGrad Steepest-Descent-like
    void useConjGradFR() { _cgmode = CGMode::CGFR; }  // reset to default
    void useConjGradPR() { _cgmode = CGMode::CGPR; } // use Polak-Ribiere CG
    void useConjGradPR0() { _cgmode = CGMode::CGPR0; } // use PR-CG with reset
    void setCGMode(CGMode cgmode) { _cgmode = cgmode; }
    CGMode getCGMode() const { return _cgmode; }

    // Setters
    void setStepSize(double stepSize) { _stepSize = stepSize; }
    void setMaxNBracketing(int maxn_bracket) { _max_n_bracketing = std::max(1, maxn_bracket); }

    // Getters
    double getStepSize() const { return _stepSize; }
    int getMaxNBracketing() const { return _max_n_bracketing; }
};
} // namespace nfm

#endif
