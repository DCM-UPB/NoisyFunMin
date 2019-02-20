#include "nfm/NoisyFunctionValue.hpp"



// --- Operations

bool NoisyFunctionValue::operator<(const NoisyFunctionValue &val)
{
    if ( _f + 2.0*_df < val._f - 2.0*val._df )
        {
            return true;
        }
    return false;
}


bool NoisyFunctionValue::operator>(const NoisyFunctionValue &val)
{
    if ( _f - 2.0*_df > val._f + 2.0*val._df )
        {
            return true;
        }
    return false;
}


bool NoisyFunctionValue::operator==(const NoisyFunctionValue &val)
{
    if ( ( (*this)<val ) || ( (*this)>val ) )
        {
            return false;
        }
    return true;
}


NoisyFunctionValue& NoisyFunctionValue::operator=(const NoisyFunctionValue &val)
{
    _ndim=val._ndim;
    for (int i=0; i<_ndim; ++i){_x[i]=val._x[i];}
    _f=val._f;
    _df=val._df;
    return *this;
}


// --- Setters

void NoisyFunctionValue::setX(const double * x)
{
    for (int i=0; i<_ndim; ++i){ _x[i] = x[i]; }
}

void NoisyFunctionValue::setF(const double &f, const double &df)
{
    _f = f;
    _df = df;
}



// --- Constructor and destructor

NoisyFunctionValue::NoisyFunctionValue(int ndim)
{
    _ndim = ndim;
    _x = new double[_ndim];
    for (int i=0; i<_ndim; ++i){_x[i]=0.;}
    _f=0.;
    _df=0.;
}


NoisyFunctionValue::~NoisyFunctionValue()
{
    delete[] _x;
    _ndim = 0;
}
