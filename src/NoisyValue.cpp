#include "nfm/NoisyValue.hpp"

#include <cmath>

namespace nfm
{

// --- Sigma Level

// default sigmaLevel
double NoisyValue::_sigmaLevel = DEFAULT_SIGMA_LEVEL;

void NoisyValue::setSigmaLevel(const double sigmaLevel)
{
    _sigmaLevel = (sigmaLevel > 0.) ? sigmaLevel : 0.;
}

// --- Set

void NoisyValue::set(const double value, const double error)
{
    val = value;
    err = error;
}

void NoisyValue::zero()
{
    val = 0.;
    err = 0.;
}

// --- Binary operations

// Compound assigment with scalar

NoisyValue NoisyValue::operator+=(const double rhs)
{
    val += rhs;
    return *this;
}

NoisyValue NoisyValue::operator-=(const double rhs)
{
    val -= rhs;
    return *this;
}

NoisyValue NoisyValue::operator*=(const double rhs)
{
    val *= rhs;
    err *= fabs(rhs);
    return *this;
}

NoisyValue NoisyValue::operator/=(const double rhs)
{
    val /= rhs;
    err /= fabs(rhs);
    return *this;
}

// Binary operations with scalar (reuse compound assigment)

NoisyValue operator+(NoisyValue lhs, const double rhs)
{
    lhs += rhs;
    return lhs;
}

NoisyValue operator-(NoisyValue lhs, const double rhs)
{
    lhs -= rhs;
    return lhs;
}

NoisyValue operator*(NoisyValue lhs, const double rhs)
{
    lhs *= rhs;
    return lhs;
}

NoisyValue operator/(NoisyValue lhs, const double rhs)
{
    lhs /= rhs;
    return lhs;
}

// Compound assigment with other noisy value

NoisyValue NoisyValue::operator+=(const NoisyValue rhs)
{
    val += rhs.val;
    err = sqrt(err*err + rhs.err*rhs.err); // standard error propagation
    return *this;
}

NoisyValue NoisyValue::operator-=(const NoisyValue rhs)
{
    val -= rhs.val;
    err = sqrt(err*err + rhs.err*rhs.err); // standard error propagation
    return *this;
}

// Binary operations with other noisy value (reuse compound assigment)

NoisyValue operator+(NoisyValue lhs, const NoisyValue rhs)
{
    lhs += rhs;
    return lhs;
}

NoisyValue operator-(NoisyValue lhs, const NoisyValue rhs)
{
    lhs -= rhs;
    return lhs;
}

// --- Comparison

bool NoisyValue::operator<(const double value) const
{
    return this->getUBound() < value;
}


bool NoisyValue::operator>(const double value) const
{
    return this->getLBound() > value;
}


bool NoisyValue::operator==(const double value) const
{
    return !(*this < value || *this > value);
}

bool NoisyValue::operator<(const NoisyValue other) const
{
    return this->getUBound() < other.getLBound();
}


bool NoisyValue::operator>(const NoisyValue other) const
{
    return this->getLBound() > other.getUBound();
}


bool NoisyValue::operator==(const NoisyValue other) const
{
    return !(*this < other || *this > other);
}

// --- Other

// Stream I/O
std::ostream &operator<<(std::ostream &os, const NoisyValue nv)
{   // write NoisyValue to ostream
    os << nv.val << " +- " << nv.err;
    return os;
}

// Minimal Distance
double NoisyValue::minDist(NoisyValue other) const
{
    return fabs(this->val - other.val) - _sigmaLevel*(this->err + other.err);
}
} // namespace nfm