#ifndef NFM_CONJGRAD_HPP
#define NFM_CONJGRAD_HPP

#include "nfm/NoisyFunMin.hpp"
#include "nfm/NoisyFunction.hpp"
#include "nfm/LineSearch.hpp"

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
// Works when target/gradient noisy is small,
// otherwise use a proper stochastic optimizer.
class ConjGrad: public NFM
{
private:
    MLMParams _mlmParams;  // line search configuration
    CGMode _cgmode; // which gradients to use

    // --- Internal methods
    void _findNextX(const std::vector<double> &dir); // do line-search
    void _writeCGDirectionToLog(const std::vector<double> &dir, const std::string &name) const;

    // --- Minimization
    void _findMin() final; // perform noisy CG minimization

public:
    explicit ConjGrad(NoisyFunctionWithGradient * targetfun, MLMParams params = defaultMLMParams());
    ~ConjGrad() final = default;

    // CG Configuration
    void useRawGrad() { _cgmode = CGMode::NOCG; }  // make ConjGrad Steepest-Descent-like
    void useConjGradFR() { _cgmode = CGMode::CGFR; }  // reset to default
    void useConjGradPR() { _cgmode = CGMode::CGPR; } // use Polak-Ribiere CG
    void useConjGradPR0() { _cgmode = CGMode::CGPR0; } // use PR-CG with reset
    void setCGMode(CGMode cgmode) { _cgmode = cgmode; }
    CGMode getCGMode() const { return _cgmode; }

    // Setters
    void setMLMParams(MLMParams params) { _mlmParams = params; }
    void setStepSize(double stepSize) { _mlmParams.stepRight = stepSize; }
    void setBackStep(double backStep) { _mlmParams.stepLeft = backStep; }
    void setMaxNBracket(int maxn_bracket) { _mlmParams.maxNBracket = maxn_bracket; }
    void setMaxNMin1D(int maxn_min1d) { _mlmParams.maxNMinimize = maxn_min1d; }

    // Getters
    MLMParams & getMLMParams() { return _mlmParams; }
    const MLMParams & getMLMParams() const { return _mlmParams; }
    double getStepSize() const { return _mlmParams.stepRight; }
    double getBackStep() const { return _mlmParams.stepLeft; }
    int getMaxNBracket() const { return _mlmParams.maxNBracket; }
    int setMaxNMin1D() const { return _mlmParams.maxNMinimize; }
};
} // namespace nfm

#endif
