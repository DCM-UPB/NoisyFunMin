#include "nfm/LineSearch.hpp"
#include "nfm/LogManager.hpp"

#include <cmath>
#include <functional>
#include <iostream>
#include <sstream>


// --- Internal Functions

// - Non-throwing checks

// check bracket for position tolerances:
// - epsx: Minimum bracket size / numerical tolerance
inline bool checkBracketXTol(const nfm::NoisyBracket &bracket, const double epsx)
{
    return (fabs(bracket.c.x - bracket.b.x) > epsx && fabs(bracket.b.x - bracket.a.x) > epsx); // simply check for the position distances
    //return fabs(bracket.c.x - bracket.a.x) > epsx*0.5*fabs(bracket.c.x + bracket.a.x) + epsx; // standard tolerance check
}

// check bracket for function value tolerances (noisy version)
// - epsf: Minimal noisy value distance between a<->b or b<->c
inline bool checkBracketFTol(const nfm::NoisyBracket &bracket, const double epsf)
{
    const double fdist = (bracket.a.f < bracket.c.f)
                         ? (fabs(bracket.a.f.value - bracket.b.f.value) - bracket.a.f.error - bracket.b.f.error)
                         : (fabs(bracket.c.f.value - bracket.b.f.value) - bracket.c.f.error - bracket.b.f.error);
    return fdist > epsf;
}

// are there neighbouring values that are equal (within error) ?
inline bool hasEquals(const nfm::NoisyBracket &bracket)
{
    return ((bracket.a.f == bracket.b.f) || (bracket.b.f == bracket.c.f));
}

// does bracket fulfill the bracketing condition a.f. > b.f < c.f
inline bool isBracketed(const nfm::NoisyBracket &bracket)
{
    return (bracket.a.f > bracket.b.f && bracket.b.f < bracket.c.f);
}

// - Throwing checks

// throw on invalid bracket X
inline void validateBracketX(const double ax, const double bx, const double cx, const std::string &callerName)
{
    if (ax >= cx) {
        throw std::invalid_argument("[" + callerName + "->validateBracketX] Bracket violates (a.x < c.x).");
    }
    if (bx >= cx || bx <= ax) {
        throw std::invalid_argument("[" + callerName + "->validateBracketX] Bracket violates (a.x < b.x < c.x).");
    }
}

// throw on invalid bracket
inline void validateBracket(const nfm::NoisyBracket &bracket, const std::string &callerName)
{
    validateBracketX(bracket.a.x, bracket.b.x, bracket.c.x, callerName);
    if (!isBracketed(bracket)) {
        throw std::invalid_argument("[" + callerName + "->validateBracket] Bracket violates (a.f > b.f < c.f).");
    }
}

// Other helper functions
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

bool findBracket(NoisyFunction &f1d, NoisyBracket &bracket /*inout*/, const int maxNIter, double epsx)
{
    // Noisy findBracket Algorithm
    // Mix of own ideas and GSL's findBracket ( gsl/min/bracketing.c )
    // Returns true when valid bracket found, else false.
    using namespace m1d_detail;

    // --- Sanity
    if (f1d.getNDim() != 1) {
        throw std::invalid_argument("[nfm::findBracket] The NoisyFunction is not 1D. Ndim=" + std::to_string(f1d.getNDim()));
    }
    bracket = sortedBracket(bracket); // ensure proper ordering
    validateBracketX(bracket.a.x, bracket.b.x, bracket.c.x, "nfm::findBracket"); // ensure valid bracket (else throw)
    epsx = std::max(0., epsx);

    // --- Initialization
    NoisyIOPair1D &a = bracket.a;
    NoisyIOPair1D &b = bracket.b;
    NoisyIOPair1D &c = bracket.c;

    // shortcut lambda
    int iter = 0; // keeps track of loop iteration count (overall)
    std::vector<double> xvec(1); // helper array to invoke noisy function
    std::function<NoisyValue(double x)> F = [&](const double x)
    {
        xvec[0] = x;
        return f1d(xvec);
    };

    // --- Bracketing

    // Pre-Processing
    writeBracketToLog("findBracket init", bracket);
    while (hasEquals(bracket)) { // we need larger interval
        // check stopping conditions
        if (!checkBracketXTol(bracket, epsx)) { return false; }
        if (isBracketed(bracket)) {
            writeBracketToLog("findBracket final", bracket);
            return true; // return with early success
        }
        if (iter++ > maxNIter) { return false; } // evaluation limit

        // scale up
        b = c;
        c.x = a.x + (b.x - a.x)/IGOLD2;
        c.f = F(c.x);

        writeBracketToLog("findBracket pre-step", bracket);
    }

    // Main Loop
    while (!hasEquals(bracket)) { // stop if we have equal function values again
        // check other stopping conditions
        if (!checkBracketXTol(bracket, epsx)) { return false; } // bracket violates tolerances
        if (isBracketed(bracket)) {
            writeBracketToLog("findBracket final", bracket);
            return true; // return with success (i.e. a.f > b.f < c.f ruled out below)
        }
        if (iter++ > maxNIter) { return false; } // evaluation limit

        // regular iteration (equals ruled out)
        if (b.f < a.f) { // -> a.f > b.f > c.f (else we would have returned successful)
            // move up
            shiftABC(a.x, b.x, c.x, (c.x - b.x)/IGOLD2 + b.x);
            shiftABC(a.f, b.f, c.f, F(c.x));
        }
        else { // -> a.f < b.f < c.f || a.f < b.f > c.f
            // contract
            c = b;
            b.x = (c.x - a.x)*IGOLD2 + a.x;
            b.f = F(b.x);
        }

        writeBracketToLog("findBracket step", bracket);
    }
    return false;
}


NoisyIOPair1D brentMin(NoisyFunction &f1d, NoisyBracket bracket, const int maxNIter, double epsx, double epsf)
{
    //
    // Adaption of GNU Scientific Libraries's Brent minimization ( gsl/min/brent.c )
    //
    using namespace m1d_detail;

    // Sanity
    if (f1d.getNDim() != 1) {
        throw std::invalid_argument("[nfm::brentMin] The NoisyFunction is not 1D. Ndim=" + std::to_string(f1d.getNDim()));
    }
    validateBracket(bracket, "nfm::brentMin"); // check for valid bracket
    epsx = std::max(0., epsx);
    epsf = std::max(0., epsf);

    // --- Initialization

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
    v.x = lb.x + IGOLD2*(ub.x - lb.x);
    v.f = F(v.x);
    w = v;

    // --- Main Brent Loop
    for (int it = 0; it < maxNIter; ++it) {
        if (!checkBracketXTol(bracket, epsx)) { break; } // bracket size too small, return early
        if (!checkBracketFTol(bracket, epsf)) { break; } // values too close, return early (noisy version)

        const double mtolb = m.x - lb.x;
        const double mtoub = ub.x - m.x;
        const double xm = 0.5*(lb.x + ub.x);
        const double tol = 1.5e-08*fabs(m.x); // tolerance for strategy choice

        std::swap(d, e); // not sure why though
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
            double t2 = 2.*tol;
            d = p/q;
            u.x = m.x + d;
            if ((u.x - lb.x) < t2 || (ub.x - u.x) < t2) {
                d = (m.x < xm) ? tol : -tol;
            }
        }
        else {
            e = (m.x < xm) ? ub.x - m.x : -(m.x - lb.x);
            d = IGOLD2*e;
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

    // Return pair with best upper error bound (out of all stored pairs)
    // This might seem inconsistent, but has proven worthy in practice.
    writeBracketToLog("brentMin final", bracket);
    NoisyIOPair1D *pairs[5]{&m, &v, &w, &ub, &lb};
    return **std::min_element(pairs, pairs + 5, [](NoisyIOPair1D * a, NoisyIOPair1D * b) { return a->f.getUBound() < b->f.getUBound(); });
}


NoisyIOPair multiLineMin(NoisyFunction &mdf, NoisyIOPair p0Pair, const std::vector<double> &dir, MLMParams params)
{
    using namespace m1d_detail;
    // Sanity
    if (static_cast<size_t>(mdf.getNDim()) != p0Pair.x.size() || p0Pair.x.size() != dir.size()) {
        throw std::invalid_argument("[nfm::multiLineMin] The passed function and positions are inconsistent in size.");
    }
    if (params.stepLeft < 0. || params.stepRight <= 0.) {
        throw std::invalid_argument("[nfm::multiLineMin] stepLeft and stepRight must be non-negative (stepRight strictly positive).");
    }
    // these should be non-zero
    params.epsx = (params.epsx > 0) ? params.epsx : STD_XTOL;
    params.epsf = (params.epsf > 0) ? params.epsf : STD_FTOL;

    // project the original multi-dim function into a one-dim function
    FunProjection1D proj1d(&mdf, p0Pair.x, dir);

    // prepare initial bracket (allow backstep via stepLeft)
    const double ax = -params.stepLeft;
    const double cx = params.stepRight;
    const double bx = ax + (cx - ax)*IGOLD2; // golden section
    NoisyBracket bracket{{ax, (fabs(ax) == 0.) ? p0Pair.f : proj1d(ax)}, // avoid recomputation
                         {bx, proj1d(bx)},
                         {cx, proj1d(cx)}};

    if (findBracket(proj1d, bracket, params.maxNBracket, params.epsx)) { // valid bracket was stored in bracket
        // now do line-minimization via brent
        NoisyIOPair1D min1D = brentMin(proj1d, bracket, params.maxNMinimize, params.epsx, params.epsf);

        if (min1D.f <= p0Pair.f) { // reject new values that are truly larger
            // return new NoisyIOPair
            p0Pair.f = min1D.f; // store the minimal f value
            proj1d.getVecFromX(min1D.x, p0Pair.x); // get the true x position
            return p0Pair;
        }
    }
    // return the old state unchanged
    return p0Pair;
}
} // namespace nfm
