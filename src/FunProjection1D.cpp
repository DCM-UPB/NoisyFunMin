#include "FunProjection1D.hpp"



void FunProjection1D::getVecFromX(const double &x, double *vec)
{
    for (int i=0; i<_originalndim; ++i)
        {
            vec[i]=_p0[i]+_direction[i]*x;
        }
}


void FunProjection1D::f(const double *x, double &f, double &df)
{
    this->getVecFromX(*x,_vec);
    _mdf->f(_vec,f,df);
}


FunProjection1D::~FunProjection1D()
{
    _ndim=0;
    _originalndim=0;
    delete [] _p0;
    delete [] _direction;
    delete [] _vec;
    _mdf=0;
}
