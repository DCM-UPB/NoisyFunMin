#include "nfm/NoisyValue.hpp"

#include <cmath>

namespace nfm
{

// --- Sigma Level

double NoisyValue::getSigmaLevel()
{
    return (_sigmaLevel > 0.) ? _sigmaLevel : DEFAULT_SIGMA_LEVEL;
}

void NoisyValue::setSigmaLevel(const double sigmaLevel)
{
    _sigmaLevel = (sigmaLevel > 0.) ? sigmaLevel : DEFAULT_SIGMA_LEVEL;
}

// --- Set

void NoisyValue::set(const double val, const double err)
{
    value = val;
    error = err;
}

// --- Binary operations

// Compound assigment with scalar

NoisyValue NoisyValue::operator+=(const double rhs)
{
    value += rhs;
    return *this;
}

NoisyValue NoisyValue::operator-=(const double rhs)
{
    value -= rhs;
    return *this;
}

NoisyValue NoisyValue::operator*=(const double rhs)
{
    value *= rhs;
    error *= rhs;
    return *this;
}

NoisyValue NoisyValue::operator/=(const double rhs)
{
    value /= rhs;
    error /= rhs;
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
    value += rhs.value;
    error = sqrt(error*error + rhs.error*rhs.error); // standard error propagation
    return *this;
}

NoisyValue NoisyValue::operator-=(const NoisyValue rhs)
{
    value -= rhs.value;
    error = sqrt(error*error + rhs.error*rhs.error); // standard error propagation
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

bool NoisyValue::operator<(const double val) const
{
    return this->getUBound() < val;
}


bool NoisyValue::operator>(const double val) const
{
    return this->getLBound() > val;
}


bool NoisyValue::operator==(const double val) const
{
    return !(*this < val || *this > val);
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

// --- Stream I/O

std::ostream &operator<<(std::ostream &os, const NoisyValue nv)
{   // write NoisyValue to ostream
    os << nv.value << " +- " << nv.error;
    return os;
}
} // namespace nfm