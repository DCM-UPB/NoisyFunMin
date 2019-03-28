#ifndef NFM_NOISYVALUE_HPP
#define NFM_NOISYVALUE_HPP

// You can define this externally to force a different default sigma level at compile time.
// For the meaning of sigma level, see below.

#include <ostream>

// You may define this externally to force a different default sigma level
// already at compile time. For the meaning of sigma level, see below.
#ifndef DEFAULT_SIGMA_LEVEL
#define DEFAULT_SIGMA_LEVEL 1.
#endif

namespace nfm
{
// A class to represent values together with their standard error,
// assuming that the underlying distribution is normal. NoisyValues
// are the output type of the NoisyFunctions in this library.
//
// It has some overloads to aid basic arithmetic with (exact) scalar
// values and overloaded addition/substraction of two NoisyValues,
// which have a normally distributed result and automatically perform
// error propagation. Furthermore, all comparisons operators are overloaded
// with methods that use a static sigmaLevel factor, which defines within
// how many standard errors a given (exact) value is considered equal to
// a NoisyValue. Comparisons of two NoisyValues take both errors into account.
//
// NOTE 1: It was made sure that NoisyValue is an aggregate type. Among many other things,
// this means you can simply use aggregate initialization: NoisyValue x{.value = 1., .error = 0.5}
// NOTE 2: Because having only two double fields makes NoisyValues very cheap to copy, they can
//         and should simply be passed by value where a const reference would be used otherwise.
struct NoisyValue
{
private:
    // Error width factor, accessed only via get/set methods.
    // This static value is initially DEFAULT_SIGMA_LEVEL.
    static double _sigmaLevel;

public:
    double value;
    double error; // standard error (sigma)

    // Static methods
    static void setSigmaLevel(double sigmaLevel = DEFAULT_SIGMA_LEVEL); // will set 0. if sigmaLevel <= 0
    static double getSigmaLevel() { return _sigmaLevel; };

    //Setters
    void set(double val, double err); // set both fields at once

    //Getters
    double getUBound() const { return value + error*_sigmaLevel; }
    double getLBound() const { return value - error*_sigmaLevel; }

    // Binary Operators (implemented based on compound assignments)

    NoisyValue operator+=(double rhs); // add assign with scalar value
    NoisyValue operator-=(double rhs); // sub assign with scalar value
    NoisyValue operator*=(double rhs); // mul assign with scalar value
    NoisyValue operator/=(double rhs); // div assign with scalar value

    friend NoisyValue operator+(NoisyValue lhs, double rhs); // add scalar (to value)
    friend NoisyValue operator-(NoisyValue lhs, double rhs); // sub scalar (to value)
    friend NoisyValue operator*(NoisyValue lhs, double rhs); // mul scalar (scale both fields)
    friend NoisyValue operator/(NoisyValue lhs, double rhs); // div scalar (scale both fields)

    NoisyValue operator+=(NoisyValue rhs); // add assign with other noisy value
    NoisyValue operator-=(NoisyValue rhs); // sub assign with other noisy value

    friend NoisyValue operator+(NoisyValue lhs, NoisyValue rhs); // add two noisy values
    friend NoisyValue operator-(NoisyValue lhs, NoisyValue rhs); // sub two noisy values

    // Comparison Operators
    bool operator<(double val) const;
    bool operator>=(double val) const { return !(*this < val); }
    bool operator>(double val) const;
    bool operator<=(double val) const { return !(*this > val); }
    bool operator==(double val) const;
    bool operator!=(double val) const { return !(*this == val); }

    bool operator<(NoisyValue other) const;
    bool operator>=(NoisyValue other) const { return !(*this < other); }
    bool operator>(NoisyValue other) const;
    bool operator<=(NoisyValue other) const { return !(*this > other); }
    bool operator==(NoisyValue other) const;
    bool operator!=(NoisyValue other) const { return !(*this == other); }

    // Stream Output
    friend std::ostream &operator<<(std::ostream &os, NoisyValue nv);
};
} // namespace nfm

#endif
