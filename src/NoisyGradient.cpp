#include "nfm/NoisyGradient.hpp"

#include <stdexcept>
#include <cmath>

namespace nfm
{

NoisyGradient::NoisyGradient(const int ndim, const bool flag_err) {
    if (ndim<=0) {
        throw std::invalid_argument("[NoisyGradient] Number of dimensions must be at least 1.");
    }
    val.assign(static_cast<size_t>(ndim), 0.);
    if (flag_err) {
        err = val;
    }
}

bool NoisyGradient::operator>(const double value) const {
    if (this->hasErr()) {
        for (size_t i=0; i<val.size(); ++i) {
            if (fabs(val[i]) - NoisyValue::getSigmaLevel()*err[i] > fabs(value)) {
                return true;
            }
        }
    }
    else {
        for (auto & gv : val) {
            if (fabs(gv) > fabs(value)) {
                return true;
            }
        }
    }
    return false;
}

void NoisyGradient::set(const int i, const NoisyValue nv) {
    val[i] = nv.val;
    err[i] = nv.err;
}
} // namespace nfm