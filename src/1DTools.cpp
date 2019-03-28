#include "nfm/1DTools.hpp"
#include "nfm/LogManager.hpp"

#include <cmath>
#include <functional>
#include <iostream>
#include <sstream>


// --- Internal Functions

// check bracket for position tolerances:
// - epsx: Minimum bracket size / numerical tolerance
inline bool checkBracketXTol(const nfm::NoisyBracket &bracket, const double epsx)
{
    if (bracket.c.x - bracket.a.x < epsx) { return false; } // simply check for the bracket size
    return bracket.c.x - bracket.a.x >= epsx*((bracket.c.x + bracket.a.x)*0.5) + epsx; // numerical tolerance check
}

// check bracket for function value tolerances (noisy version)
// - epsf: Minimal noisy value distance between a<->b or b<->c
inline bool checkBracketFTol(const nfm::NoisyBracket &bracket, const double epsf)
{
    const double fdist = (bracket.a.f < bracket.c.f)
                            ? (fabs(bracket.a.f.value - bracket.b.f.value) - bracket.a.f.error - bracket.b.f.error)
                            : (fabs(bracket.c.f.value - bracket.b.f.value) - bracket.c.f.error - bracket.b.f.error);
    return fdist >= epsf;
}

// throw on invalid bracket X
inline void checkBracketX(const double ax, const double bx, const double cx, const std::string &callerName)
{
    if (ax >= cx) {
        throw std::invalid_argument("[" + callerName + "->checkBracketX] Bracket violates (a.x < c.x).");
    }
    if (bx >= cx || bx <= ax) {
        throw std::invalid_argument("[" + callerName + "->checkBracketX] Bracket violates (a.x < b.x < c.x).");
    }
}

// throw on invalid bracket
inline void checkBracket(const nfm::NoisyBracket &bracket, const std::string &callerName)
{
    checkBracketX(bracket.a.x, bracket.b.x, bracket.c.x, callerName);
    if (bracket.b.f >= bracket.a.f || bracket.b.f >= bracket.c.f) {
        throw std::invalid_argument("[" + callerName + "->checkBracket] Bracket violates (a.f > b.f < c.x).");
    }
}

// other helper functions
template <class T>
inline void shiftABC(T &a, T &b, T &c, const T d)
{
    a = b;
    b = c;
    c = d;
}

inline nfm::NoisyBracket sortedBracket(nfm::NoisyBracket bracket) // we take value and return by value
{   // make sure bracket x's are in ascending order (assuming b is bracketed)
    if (bracket.a.x > bracket.c.x) { std::swap(bracket.a, bracket.c); }
    return bracket;
}

// --- Public Functions

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

bool findBracket(NoisyFunction &f1d, NoisyBracket &bracket /*inout*/, double epsx)
{
    //
    // Adaption of GNU Scientific Libraries's findBracket ( gsl/min/bracketing.c )
    //

    // Sanity
    if (f1d.getNDim() != 1) {
        throw std::invalid_argument("[nfm::findBracket] The NoisyFunction is not 1D. Ndim=" + std::to_string(f1d.getNDim()));
    }
    bracket = sortedBracket(bracket); // ensure proper ordering
    checkBracketX(bracket.a.x, bracket.b.x, bracket.c.x, "nfm::findBracket"); // ensure valid bracket
    epsx = std::max(0., epsx);

    // Initialization
    constexpr double GOLDEN = 0.3819660; // stretch factor for successive steps
    constexpr int NEVAL_LIMIT = m1d_default::MAX_NEVAL; // how many function evaluations to allow until throwing

    NoisyIOPair1D &a = bracket.a;
    NoisyIOPair1D &b = bracket.b;
    NoisyIOPair1D &c = bracket.c;

    // shortcut lambda
    std::vector<double> xvec(1); // helper array to invoke noisy function
    std::function<NoisyValue(double x)> F = [&](const double x)
    {
        xvec[0] = x;
        return f1d(xvec);
    };

    // Adaption of GSL's findBracket algorithm

    // initial step
    int ieval = 1; // keeps track of number of function evaluations ( we count the initial one already)
    if (c.f.value >= a.f.value) { // no reason to involve errors in this decision
        b.x = (c.x - a.x)*GOLDEN + a.x;
        b.f = F(b.x);
    }
    else {
        b = c;
        c.x = (b.x - a.x)/GOLDEN + a.x;
        c.f = F(c.x);
    }
    writeBracketToLog("findBracket init", bracket);

    // Main findBracket loop. If bracket OK, returns true.
    // If false is returned, no valid bracket could be found.
    while (true) {
        if (a.f > b.f && b.f < c.f) { // return with success
            writeBracketToLog("findBracket final", bracket);
            return true;
        }
        else { // check tolerances, if fail return without success
            if (ieval++ > NEVAL_LIMIT) { return false; } // evaluation limit
            if (!checkBracketXTol(bracket, epsx)) { return false; } // bracket violates tolerances
        }

        // continue with the iteration
        if (b.f < a.f) {
            if (b.f > c.f) {
                shiftABC(a.x, b.x, c.x, (b.x - a.x)/GOLDEN + a.x);
                shiftABC(a.f, b.f, c.f, F(c.x));
            }
            else { // b.f == c.f
                c = b;
                b.x = (c.x - a.x)*GOLDEN + a.x;
                b.f = F(b.x);
            }
        }
        else { // b.f >= a.f
            c = b;
            b.x = (c.x - a.x)*GOLDEN + a.x;
            b.f = F(b.x);
        }

        writeBracketToLog("findBracket step", bracket);
    }
}


NoisyIOPair1D brentMin(NoisyFunction &f1d, NoisyBracket bracket, double epsx, double epsf)
{
    //
    // Adaption of GNU Scientific Libraries's Brent minimization ( gsl/min/brent.c )
    //

    // Sanity
    if (f1d.getNDim() != 1) {
        throw std::invalid_argument("[nfm::brentMin] The NoisyFunction is not 1D. Ndim=" + std::to_string(f1d.getNDim()));
    }
    checkBracket(bracket, "nfm::brentMin"); // check for valid bracket
    epsx = std::max(0., epsx);
    epsf = std::max(0., epsf);

    // Initialization
    constexpr int MAX_NITERATIONS = m1d_default::MAX_NEVAL;
    constexpr double GOLDEN = 0.3819660;

    // we reuse the bracket
    NoisyIOPair1D &lb = bracket.a; // lower bound
    NoisyIOPair1D &m = bracket.b;
    NoisyIOPair1D &ub = bracket.c; // upper bound

    // shortcut lambda
    std::vector<double> xvec(1); // helper array to invoke noisy function
    std::function<NoisyValue(double x)> F = [&](const double x)
    {
        xvec[0] = x;
        return f1d(xvec);
    };

    // initialize helpers
    double d = 0., e = 0.;
    NoisyIOPair1D v{}, w{};
    v.x = lb.x + GOLDEN*(ub.x - lb.x);
    v.f = F(v.x);
    w = v;

    // main brent loop
    for (int it = 0; it < MAX_NITERATIONS; ++it) {
        if (!checkBracketXTol(bracket, epsx)) { break; } // bracket size too small, return early
        if (!checkBracketFTol(bracket, epsf)) { break; } // values too close, return early (noisy version)

        const double mtolb = m.x - lb.x;
        const double mtoub = ub.x - m.x;
        const double xm = 0.5*(lb.x + ub.x);
        const double tol = 1.5e-08*fabs(m.x); // numeric tolerance

        std::swap(d, e); // not sure why this is in the original
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
            d = GOLDEN*e;
        }

        if (fabs(d) >= tol) {
            u.x = m.x + d;
        }
        else {
            u.x = m.x + ((d > 0) ? tol : -tol);
        }

        // here we evaluate the function
        u.f = F(u.x);

        // check continue conditions
        if (u.f <= m.f) {
            if (u.x < m.x) { ub = m; }
            else { lb = m; }

            v = w;
            w = m;
            m = u;
            writeBracketToLog("brentMin step", bracket);
            continue;
        }

        if (u.x < m.x) { lb = u; }
        else { ub = u; }

        if (u.f <= w.f || w.x == m.x) {
            v = w;
            w = u;
            writeBracketToLog("brentMin step", bracket);
            continue;
        }
        if (u.f <= v.f || v.x == m.x || v.x == w.x) {
            v = u;
            writeBracketToLog("brentMin step", bracket);
            continue;
        }

        writeBracketToLog("brentMin step", bracket);
    }

    // return best result
    return m;
}

NoisyIOPair multiLineMin(NoisyFunction &mdf, NoisyIOPair p0Pair, const std::vector<double> &dir,
                         double stepLeft, double stepRight, double epsx, double epsf)
{
    // Sanity
    if (mdf.getNDim() != p0Pair.x.size() || mdf.getNDim() != dir.size()) {
        throw std::invalid_argument("[nfm::multiLineMin] The passed function and positions are inconsistent in size.");
    }
    if (stepLeft <= 0. || stepRight <= 0.) {
        throw std::invalid_argument("[nfm::multiLineMin] stepLeft and stepRight must be positive numbers.");
    }

    // project the original multi-dim function into a one-dim function
    FunProjection1D proj1d(&mdf, p0Pair.x, dir);

    // prepare bracket boundaries (allow backstep via stepLeft)
    NoisyBracket bracket{{-stepLeft, proj1d(-stepLeft)},
                         {0.,              p0Pair.f},
                         {stepRight,       proj1d(stepRight)}};

    // find initial bracket and then the minimum in the bracket
    if (findBracket(proj1d, bracket, epsx)) { // valid bracket was stored in bracket
        // line-minmization via brent
        NoisyIOPair1D min1D = brentMin(proj1d, bracket, epsx, epsf);

        // return NoisyIOPair
        p0Pair.f = min1D.f; // store the minimal f value
        proj1d.getVecFromX(min1D.x, p0Pair.x); // get the true x position
    }

    // return the new state (unchanged if bracketing failed)
    return p0Pair;
}
} // namespace nfm
