#ifndef NFM_FUNPROJECTION1D_HPP
#define NFM_FUNPROJECTION1D_HPP

#include "nfm/NoisyFunction.hpp"

#include <algorithm>

namespace nfm
{

class FunProjection1D final: public NoisyFunction
{
private:
    NoisyFunction * const _mdf;  //multidimensional function that must be projected
    const std::vector<double> _p0;   //starting point
    const std::vector<double> _dir;   //direction
    std::vector<double> _vec;  //vector used internally

public:
    FunProjection1D(NoisyFunction * mdf, std::vector<double> p0, std::vector<double> dir);
    ~FunProjection1D() final = default;

    // calculate true vector from scalar projection coordinate
    void getVecFromX(double x, std::vector<double> &vec /*out*/);

    //projected one-dimensional function
    NoisyValue f(const std::vector<double> &x) final;

    //projected one-dimensional function (using single double)
    NoisyValue f(double x);
    NoisyValue operator()(double x) { return this->f(x); }
};
} // namespace nfm

#endif
