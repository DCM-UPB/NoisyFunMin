#include "nfm/FunProjection1D.hpp"

#include <iostream>

namespace nfm
{

FunProjection1D::FunProjection1D(NoisyFunction * mdf, std::vector<double> p0, std::vector<double> dir):
        NoisyFunction(1), _mdf(mdf), _p0(std::move(p0)), _dir(std::move(dir))
{
    if (_p0.size() != static_cast<size_t>(_mdf->getNDim())) {
        throw std::invalid_argument("[FunProjection1D] Size of the initial position vector is not equal to NoisyFunction dimension.");
    }
    if (_dir.size() != static_cast<size_t>(_mdf->getNDim())) {
        throw std::invalid_argument("[FunProjection1D] Size of the direction vector is not equal to NoisyFunction dimension.");
    }
    _vec.assign(_p0.size(), 0.);
}

void FunProjection1D::getVecFromX(const double x, std::vector<double> &vec)
{
    std::cout << "getVecFromX: " << x << std::endl;
    for (int i=0; i<_mdf->getNDim(); ++i) {
        std::cout << "i " << i << " vec[i] " << _vec[i] << " p0[i] " << _p0[i] << " dir[i] " << _dir[i] << std::endl;
        _vec[i] = _p0[i] + x*_dir[i];
        std::cout << "-> vec[i] = " << _vec[i] << std::endl;
    }
}

NoisyValue FunProjection1D::f(const std::vector<double> &x) {
    this->getVecFromX(x[0], _vec);
    return _mdf->f(_vec);
}
} // namespace nfm