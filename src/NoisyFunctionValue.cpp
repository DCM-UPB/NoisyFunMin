#include "nfm/NoisyFunctionValue.hpp"

#include <algorithm>
#include <stdexcept>

// --- Operations

bool NoisyFunctionValue::operator<(const NoisyFunctionValue &val)
{
    return _f + 2.0*_df < val._f - 2.0*val._df;
}


bool NoisyFunctionValue::operator>(const NoisyFunctionValue &val)
{
    return _f - 2.0*_df > val._f + 2.0*val._df;
}


bool NoisyFunctionValue::operator==(const NoisyFunctionValue &val)
{
    return !(( (*this)<val ) || ( (*this)>val ));
}


NoisyFunctionValue& NoisyFunctionValue::operator=(const NoisyFunctionValue &val)
{   // we actually copy here instead of assign... but will do for the moment
    if (_ndim != val._ndim) {
        throw std::runtime_error("[NoisyFunctionValue] For copying contents, we need both noisy value's ndim to be identical.");
    }
    std::copy(val._x, val._x+_ndim, _x);
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
    std::fill(_x, _x+_ndim, 0.);
    _f=0.;
    _df=0.;
}


NoisyFunctionValue::~NoisyFunctionValue()
{
    delete[] _x;
}
