#ifndef NFM_FUNPROJECTION1D_HPP
#define NFM_FUNPROJECTION1D_HPP

#include "nfm/NoisyFunction.hpp"

#include <algorithm>
#include <stdexcept>

namespace nfm
{

class FunProjection1D: public NoisyFunction
{
protected:
    NoisyFunction * const _mdf;  //multidimensional function that must be projected
    const std::vector<double> _p0;   //starting point
    const std::vector<double> _direction;   //direction
    std::vector<double> _vec;  //vector used internally

public:
    FunProjection1D(NoisyFunction * mdf, std::vector<double> p0, std::vector<double> direction);
    ~FunProjection1D() final = default;

    // calculate true vector from scalar projection coordinate
    void getVecFromX(double x, std::vector<double> &vec /*out*/);

    //projected one-dimensional function
    NoisyValue f(const std::vector<double> &x) final;
};
} // namespace nfm

#endif
