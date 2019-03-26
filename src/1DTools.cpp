#include "nfm/1DTools.hpp"

#include "nfm/LogManager.hpp"

#include <sstream>
#include <cmath>
#include <iostream>
#include <functional>


// --- Internal helpers

// called if maximum number of evaluations is reached
inline void abortFindBracket()
{
    throw std::runtime_error("NoisyFunctionMin Error! Bracketing is taking way too long...");
}

// used to evaluate F1D and increase/check evaluation count
inline nfm::NoisyValue evalF1D(const double x, std::vector<double> &xv, nfm::NoisyFunction &f1d, int &ieval, const int neval_max)
{
    if (++ieval > neval_max) { abortFindBracket(); }
    xv[0] = x;
    return f1d(xv);
}

// helper functions
template <class T>
inline void shiftABC(T &a, T &b, T &c, const T d)
{
    a = b;
    b = c;
    c = d;
}

inline nfm::NoisyBracket sortBracket(nfm::NoisyBracket bracket) // we take value and return by value
{   // make sure bracket x's are in ascending order (assuming b is bracketed)
    if (bracket.a.x > bracket.c.x) { std::swap(bracket.a, bracket.c); }
    return bracket;
}

namespace nfm
{

void writeBracketToLog(const std::string &key, const NoisyBracket &bracket)
{
    std::stringstream s;
    s << key << ":    " <<
      bracket.a.x << " -> " << bracket.a.f << "    " <<
      bracket.b.x << " -> " << bracket.b.f << "    " <<
      bracket.c.x << " -> " << bracket.c.f;
    s << std::flush;
    LogManager::logString(s.str(), LogLevel::VERBOSE);
}

NoisyBracket findBracket(NoisyFunction &f1d, const double initX1, const double initX2)
{
    if (f1d.getNDim() != 1) {
        throw std::invalid_argument("[nfm::findBracket] The NoisyFunction is not 1D. Ndim=" + std::to_string(f1d.getNDim()));
    }

    constexpr double STRETCH_FACTOR = 1.618; // stretch factor for successive steps
    constexpr double STRETCH_LIMIT = 100.; // maximal stretch factor
    constexpr double MIN_DOUBLE = 1.e-20; // prevent malicious division
    constexpr int NEVAL_LIMIT = 100; // how many function evaluations to allow until throwing

    NoisyBracket bracket{};
    NoisyIOPair1D &a = bracket.a;
    NoisyIOPair1D &b = bracket.b;
    NoisyIOPair1D &c = bracket.c;

    NoisyIOPair1D h{}; // helper variable
    std::vector<double> xvec(1); // helper array to invoke noisy function

    int ieval = 0; // keeps track of number of function evaluations
    // shortcut lambda
    std::function<NoisyValue(double x)> F = [&](const double x) { return evalF1D(x, xvec, f1d, ieval, NEVAL_LIMIT); };

    a.x = initX1;
    a.f = F(a.x);
    b.x = initX2;
    b.f = F(b.x);
    if (b.f > a.f) { std::swap(a, b); } // make sure that downhill is in direction of b

    c.x = b.x + STRETCH_FACTOR*(b.x - a.x); // start with first guess
    c.f = F(c.x);
    writeBracketToLog("findBracket init", bracket);

    while (b.f > c.f) { // loop until bracket found
        // compute parabolic extrapolation from bracket
        const double r = (b.x - a.x)*(b.f.value - c.f.value);
        const double q = (b.x - c.x)*(b.f.value - a.f.value);
        const double qmr = q - r; // used for sign
        const double qmr2 = 2.*std::max(fabs(qmr), MIN_DOUBLE); // prevent bad division
        h.x = b.x - ((b.x - c.x)*q - (b.x - a.x)*r)/(qmr >= 0. ? qmr2 : -qmr2);
        const double hlim = b.x + STRETCH_LIMIT*(c.x - b.x); // limit on hx

        // check all cases
        if ((b.x - h.x)*(h.x - c.x) > 0.) { // h is between b and c
            h.f = F(h.x);
            if (h.f < c.f) { // minimum between b and c
                a = b;
                b = h;
                writeBracketToLog("findBracket final", bracket);
                return sortBracket(bracket); // return ascending-x bracket
            }
            else if (h.f > b.f) {// minimum between a and h
                c = h;
                writeBracketToLog("findBracket final", bracket);
                return sortBracket(bracket); // return ascending-x bracket
            }
            h.x = c.x + STRETCH_FACTOR*(c.x - b.x); // parabolic fit didn't work, use default stretch
            h.f = F(h.x);
        }
        else if ((c.x - h.x)*(h.x - hlim) > 0.) { // fit is between c and hlim
            h.f = F(h.x);
            if (h.f < c.f) {
                shiftABC(b.x, c.x, h.x, h.x + STRETCH_FACTOR*(h.x - c.x));
                shiftABC(b.f, c.f, h.f, F(h.x));
            }
        }
        else if ((h.x - hlim)*(hlim - c.x) >= 0.) { // Limit h to max value
            h.x = hlim;
            h.f = F(h.x);
        }
        else { // reject h, use default stretch
            h.x = c.x + STRETCH_FACTOR*(c.x - b.x);
            h.f = F(h.x);
        }
        // prepare next iteration and continue
        shiftABC(a, b, c, h);
        writeBracketToLog("findBracket step", bracket);
    }

    writeBracketToLog("findBracket final", bracket);
    return sortBracket(bracket); // return ascending-x bracket
}

NoisyIOPair1D brentMinimization(NoisyFunction &f1d, const double eps, NoisyBracket bracket)
{
    //
    // Adaption of GNU Scientific Libraries's Brent Minimization for NoisyValues
    //
    if (f1d.getNDim() != 1) {
        throw std::invalid_argument("[nfm::brentMinimization] The NoisyFunction is not 1D. Ndim=" + std::to_string(f1d.getNDim()));
    }
    if (bracket.b.f > bracket.a.f) { throw std::invalid_argument("[brentMinimization]: Initial bracket contains fb>fa"); }
    if (bracket.b.f > bracket.c.f) { throw std::invalid_argument("[brentMinimization]: Initial bracket contains fb>fc"); }
    if (bracket.a.x > bracket.c.x) { throw std::invalid_argument("[brentMinimization]: Initial bracket contains ax>cx"); }

    // shortcut lambda
    std::vector<double> xvec(1); // helper array to invoke noisy function
    std::function<NoisyValue(double x)> F = [&](const double x)
    {
        xvec[0] = x;
        return f1d(xvec);
    };

    const int MAX_NITERATIONS = 100;
    const double GOLDEN_MEAN = 0.382;

    // we reuse the bracket
    NoisyIOPair1D &lb = bracket.a; // lower bound
    NoisyIOPair1D &m = bracket.b;
    NoisyIOPair1D &ub = bracket.c; // upper bound

    // initialize helpers
    double d = 0., e = 0.;
    NoisyIOPair1D v{}, w{};
    v.x = lb.x + GOLDEN_MEAN*(ub.x - lb.x);
    v.f = F(v.x);
    w = v;

    // main brent loop
    for (int it = 0; it < MAX_NITERATIONS; ++it) {
        const double present_eps = (lb.f > ub.f)
                                   ? (lb.f.value - m.f.value - lb.f.error - m.f.error)
                                   : (ub.f.value - m.f.value - ub.f.error - m.f.error);
        if (present_eps > eps) { break; } // terminate early

        const double mtolb = m.x - lb.x;
        const double mtoub = ub.x - m.x;
        const double xm = 0.5*(lb.x + ub.x);
        const double tol = 1.5e-08*fabs(m.x); // numeric tolerance

        std::swap(d, e);
        NoisyIOPair1D u{};

        double p = 0.;
        double q = 0.;
        double r = 0.;

        // fit parabola (magic code from GSL)
        if (fabs(e) > tol) {
            r = (m.x - w.x)*(m.f.value - v.f.value);
            q = (m.x - v.x)*(m.f.value - w.f.value);
            p = (m.x - v.x)*q - (m.x - w.x)*r;
            q = 2.*(q - r);

            if (q > 0.) {
                p = -p;
            }
            else {
                q = -q;
            }

            r = e;
            e = d;
        }

        if (fabs(p) < fabs(0.5*q*r) && p < q*mtolb && p < q*mtoub) {
            double t2 = 2*tol;

            d = p/q;
            u.x = m.x + d;

            if ((u.x - lb.x) < t2 || (ub.x - u.x) < t2) {
                d = (m.x < xm) ? tol : -tol;
            }
        }
        else {
            e = (m.x < xm) ? ub.x - m.x : -(m.x - lb.x);
            d = GOLDEN_MEAN*e;
        }

        if (fabs(d) >= tol) {
            u.x = m.x + d;
        }
        else {
            u.x = m.x + ((d > 0) ? tol : -tol);
        }

        // here the function get's evaluated
        u.f = F(u.x);

        // check continue conditions
        if (u.f <= m.f) {
            if (u.x < m.x) { ub = m; }
            else { lb = m; }

            v = w;
            w = m;
            m = u;
            continue;
        }
        else {
            if (u.x < m.x) { lb = u; }
            else { ub = u; }

            if (u.f <= w.f || w.x == m.x) {
                v = w;
                w = u;
                continue;
            }
            else if (u.f <= v.f || v.x == m.x || v.x == w.x) {
                v = u;
                continue;
            }
        }
    }

    // return best result
    return m;
}

NoisyIOPair multiLineMinimization(NoisyFunction &mdf, const std::vector<double> &p0, const std::vector<double> &dir, double eps, double initX1, double initX2)
{
    // project the original multi-dim function into a one-dim function
    FunProjection1D proj1d(&mdf, p0, dir);

    // find initial bracket and then the minimum in the bracket
    NoisyIOPair1D min1D = brentMinimization(proj1d, eps, findBracket(proj1d, initX1, initX2));

    // return NoisyIOPair
    NoisyIOPair minND;
    minND.f = min1D.f; // store the minimal f value
    proj1d.getVecFromX(min1D.x, minND.x); // get the true x position
    return minND;
}
} // namespace nfm