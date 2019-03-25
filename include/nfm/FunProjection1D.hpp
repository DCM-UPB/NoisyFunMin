#ifndef NFM_FUNPROJECTION1D_HPP
#define NFM_FUNPROJECTION1D_HPP

#include "nfm/NoisyFunction.hpp"

#include <algorithm>

namespace nfm
{

class FunProjection1D: public NoisyFunction
{
protected:
    double * _p0;   //starting point
    double * _direction;   //direction
    double * _vec;  //vector used internally
    NoisyFunction * _mdf;  //multidimensional function that must be projected
    int _originalndim;  //dimension of the multidimensional function

public:
    FunProjection1D(const int & originalndim, const double * p0, const double * direction, NoisyFunction * mdf): NoisyFunction(1)
    {
        _originalndim=originalndim;
        _p0=new double[_originalndim];
        _direction=new double[_originalndim];
        _vec=new double[_originalndim];

        std::copy(p0, p0+_originalndim, _p0);
        std::copy(direction, direction+_originalndim, _direction);
        std::fill(_vec, _vec+_originalndim, 0.);

        _mdf=mdf;
    }
    ~FunProjection1D() override;

    void f(const double *x, double &f, double &df) override;  //projected one-dimensional function

    void getVecFromX(const double &x, double *vec);
};
} // namespace nfm

#endif
