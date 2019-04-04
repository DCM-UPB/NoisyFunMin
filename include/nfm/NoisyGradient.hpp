#ifndef NFM_NOISYGRADIENT_HPP
#define NFM_NOISYGRADIENT_HPP

#include "nfm/NoisyValue.hpp"

#include <vector>

namespace nfm
{

// A class to store gradients together with their error.
// Mimics NoisyValue, but has only some essential methods.
// NOTE: The contained vectors are publicly available,
//       but please never change their size (directly).
struct NoisyGradient
// Used to store noisy gradients
{
    std::vector<double> val; //
    std::vector<double> err;

    explicit NoisyGradient(int ndim, bool flag_err);

    // get the dimensions
    int getNDim() const { return static_cast<int>(val.size()); }
    size_t size() const { return val.size(); }
    bool empty() const { return val.empty(); }
    bool hasErr() const { return !err.empty(); }

    // set/get elements in NoisyValue view
    void set(int i, NoisyValue nv);
    NoisyValue get(int i) const { return {val[i], err[i]}; }
    const NoisyValue operator[](size_t i) const { return {val[i], err[i]}; }

    // Specialized Comparison:
    // Does at least one element gi fulfill: abs(gi_val) - gi_err > abs(val) ?
    bool operator>(double value) const;
    bool operator<(double value) const { return !(*this > value); }
};
} // namespace nfm

#endif
