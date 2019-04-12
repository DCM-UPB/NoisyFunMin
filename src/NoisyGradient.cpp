#include "nfm/NoisyGradient.hpp"

#include <cmath>

namespace nfm
{

NoisyGradient::NoisyGradient(const int ndim)
{
    if (ndim <= 0) {
        throw std::invalid_argument("[NoisyGradient] Number of dimensions must be at least 1.");
    }
    val.assign(static_cast<size_t>(ndim), 0.);
    err.assign(val.size(), 0.);
}

void NoisyGradient::zero()
{
    std::fill(val.begin(), val.end(), 0.);
    std::fill(err.begin(), err.end(), 0.);
}

void NoisyGradient::set(const NoisyValue nv)
{
    std::fill(val.begin(), val.end(), nv.val);
    std::fill(err.begin(), err.end(), nv.err);
}

void NoisyGradient::set(const int i, const NoisyValue nv)
{
    val[i] = nv.val;
    err[i] = nv.err;
}

bool NoisyGradient::operator>(const double value) const
{
    for (size_t i = 0; i < val.size(); ++i) {
        if (fabs(val[i]) - NoisyValue::getSigmaLevel()*err[i] > fabs(value)) {
            return true;
        }
    }
    return false;
}
} // namespace nfm