#include "nfm/FunProjection1D.hpp"

namespace nfm
{

FunProjection1D::FunProjection1D(NoisyFunction * mdf, std::vector<double> p0, std::vector<double> direction):
        NoisyFunction(1), _mdf(mdf), _p0(std::move(p0)), _direction(std::move(direction))
{
    if (_p0.size() != _mdf->getNDim()) {
        throw std::invalid_argument("[FunProjection1D] Size of the initial position vector is not equal to NoisyFunction dimension.");
    }
    if (_direction.size() != _mdf->getNDim()) {
        throw std::invalid_argument("[FunProjection1D] Size of the direction vector is not equal to NoisyFunction dimension.");
    }
    _vec.reserve(p0.size());
}

void FunProjection1D::getVecFromX(const double x, std::vector<double> &vec)
{
    for (int i=0; i<_mdf->getNDim(); ++i) {
        _vec[i] = _p0[i] + x*_direction[i];
    }
}

NoisyValue FunProjection1D::f(const std::vector<double> &x) {
    this->getVecFromX(x[0], _vec);
    return _mdf->f(_vec);
}
} // namespace nfm