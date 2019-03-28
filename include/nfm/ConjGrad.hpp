#ifndef NFM_CONJGRAD_HPP
#define NFM_CONJGRAD_HPP

#include "nfm/NoisyFunMin.hpp"
#include "nfm/NoisyFunction.hpp"

namespace nfm
{

class ConjGrad: public NFM
{
private:
    bool _use_conjgrad;

protected:
    // --- Internal methods
    void _findNextX(const std::vector<double> &dir); // do line-search
    void _writeCGDirectionToLog(const std::vector<double> &dir, const std::string &name) const;

    // --- Minimization
    void _findMin() final; // perform noisy CG minimization

public:
    explicit ConjGrad(NoisyFunctionWithGradient * targetfun, const int &max_n_const_values = 20):
            NFM(targetfun, max_n_const_values), _use_conjgrad(true)
    {}
    ~ConjGrad() final = default;


    // Configuration
    void useSimpleGradient() { _use_conjgrad = false; }  // make ConjGrad a Steepest Descent
    void useConjugateGradient() { _use_conjgrad = true; }  // reset to normal
};
} // namespace nfm

#endif
