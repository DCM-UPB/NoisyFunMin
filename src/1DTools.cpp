#include "nfm/1DTools.hpp"
#include "nfm/LogManager.hpp"

#include <cmath>
#include <functional>
#include <iostream>
#include <sstream>


// --- Internal Functions

// Line Minimization Accuracy (noisy version)
inline double getBracketEpsF(const nfm::NoisyBracket &bracket)
{
    return (bracket.a.f < bracket.c.f)
           ? (fabs(bracket.a.f.value - bracket.b.f.value) - bracket.a.f.error - bracket.b.f.error)
           : (fabs(bracket.c.f.value - bracket.b.f.value) - bracket.c.f.error - bracket.b.f.error);
}

// check for valid bracket X
inline void checkBracketX(const double ax, const double bx, const double cx, const std::string &callerName)
{
    if (ax >= cx) {
        throw std::invalid_argument("[" + callerName + "->checkBracketX] Bracket violates (a.x < c.x).");
    }
    if (bx >= cx || bx <= ax) {
        throw std::invalid_argument("[" + callerName + "->checkBracketX] Bracket violates (a.x < b.x < c.x).");
    }
}

// check for valid bracket
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

bool findBracket(NoisyFunction &f1d, NoisyBracket &bracket /*inout*/, const double epsx)
{
    //
    // GNU Scientific Libraries's findBracket ( gsl/min/bracketing.c ), adapted for NoisyValues
    //
    if (f1d.getNDim() != 1) {
        throw std::invalid_argument("[nfm::findBracket] The NoisyFunction is not 1D. Ndim=" + std::to_string(f1d.getNDim()));
    }

    constexpr double GOLDEN = 0.3819660; // stretch factor for successive steps
    constexpr int NEVAL_LIMIT = 100; // how many function evaluations to allow until throwing

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

    // validate input
    bracket = sortedBracket(bracket); // ensure proper ordering
    checkBracketX(bracket.a.x, bracket.b.x, bracket.c.x, "nfm::findBracket"); // ensure valid bracket

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
        if (b.f < a.f) {
            if (b.f < c.f) {
                writeBracketToLog("findBracket final", bracket);
                return true;
            }
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
        if (++ieval > NEVAL_LIMIT) { return false; }

        if ((c.x - a.x) < epsx*((c.x + a.x)*0.5) + epsx) { // bracket too small, return without success
            return false;
        }
    }

/*  // Alternative approach (doesn't work too well)
    constexpr double STRETCH_FACTOR = 1.62;
    constexpr double STRETCH_LIMIT = 100.;
    constexpr double MIN_DOUBLE = 1.e-20;

    NoisyIOPair1D h {};

    b = c;
    if (b.f > a.f) { std::swap(a, b); } // make sure that downhill is in direction of b

    c.x = b.x + STRETCH_FACTOR*(b.x - a.x); // start with first guess
    c.f = F(c.x);
    writeBracketToLog("findBracket init", bracket);

    while (b.f >= c.f || b.f >= a.f) { // loop until bracket found
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
                return sortedBracket(bracket); // return ascending-x bracket
            }
            else if (h.f > b.f) {// minimum between a and h
                c = h;
                writeBracketToLog("findBracket final", bracket);
                return sortedBracket(bracket); // return ascending-x bracket
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
    return sortedBracket(bracket); // return ascending-x bracket
*/
}


NoisyIOPair1D brentMinimization(NoisyFunction &f1d, NoisyBracket bracket, const double epsx, const double epsf)
{
    //
    // GNU Scientific Libraries's Brent minimization ( gsl/min/brent.c ), adapted for NoisyValues
    //
    if (f1d.getNDim() != 1) {
        throw std::invalid_argument("[nfm::brentMinimization] The NoisyFunction is not 1D. Ndim=" + std::to_string(f1d.getNDim()));
    }
    checkBracket(bracket, "nfm::brentMinimization"); // check for valid bracket

    constexpr int MAX_NITERATIONS = 100;
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
        if ((ub.x - lb.x) < epsx*((ub.x + lb.x)*0.5) + epsx) { break; }// bracket too small, return without success
        if (getBracketEpsF(bracket) < epsf) { break; } // terminate on tolerance check (noisy version)

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

NoisyIOPair multiLineMinimization(NoisyFunction &mdf, NoisyIOPair p0Pair, const std::vector<double> &dir,
                                  const double initWidth, const double epsx, const double epsf)
{
    // project the original multi-dim function into a one-dim function
    FunProjection1D proj1d(&mdf, p0Pair.x, dir);

    // prepare bracket boundaries (allow a bit of backstep)
    NoisyBracket bracket{{-0.25*initWidth, proj1d(-initWidth)},
                         {0.,              p0Pair.f},
                         {initWidth,       proj1d(initWidth)}};

    // find initial bracket and then the minimum in the bracket
    if (findBracket(proj1d, bracket, epsx)) { // valid bracket was stored in bracket
        // line-minmization via brent
        NoisyIOPair1D min1D = brentMinimization(proj1d, bracket, epsx, epsf);

        // return NoisyIOPair
        p0Pair.f = min1D.f; // store the minimal f value
        proj1d.getVecFromX(min1D.x, p0Pair.x); // get the true x position
    }

    // return the new state (unchanged if bracketing failed)
    return p0Pair;
}
} // namespace nfm
