#include "nfm/1DTools.hpp"

#include "nfm/LogManager.hpp"

#include <sstream>
#include <cmath>
#include <iostream>
#include <functional>
#include <nfm/1DTools.hpp>


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
inline void shift2(T &a, T &b, const T c)
{
    a = b;
    b = c;
}

template <class T>
inline void shift3(T &a, T &b, T &c, const T d)
{
    a = b;
    b = c;
    c = d;
}

template <class T>
inline void moveBracket(T &a, T &b, T &c, const T d, const T e, const T f)
{
    a = d;
    b = e;
    c = f;
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

    constexpr double STRETCH_FACTOR = 1.5; // stretch factor for successive steps
    constexpr double STRETCH_LIMIT = 100.; // maximal stretch factor
    constexpr double MIN_DOUBLE = 1.e-20; // prevent malicious division
    constexpr int NEVAL_LIMIT = 100; // how many function evaluations to allow until throwing

    NoisyBracket bracket{};
    NoisyIOPair1D &a = bracket.a;
    NoisyIOPair1D &b = bracket.b;
    NoisyIOPair1D &c = bracket.c;

    NoisyIOPair1D h{}; // helper variable
    std::vector<double> xvec(1); // helper array to invoke noisy function

    // shortcut lambda
    int ieval = 0; // keeps track of number of function evaluations
    std::function<NoisyValue(double x)> fun = [&](const double x) { return evalF1D(x, xvec, f1d, ieval, NEVAL_LIMIT); };

    a.x = initX1;
    a.f = fun(a.x);
    b.x = initX2;
    b.f = fun(b.x);
    if (b.f > a.f) { std::swap(a, b); } // make sure that downhill is in direction of b

    c.x = b.x + STRETCH_FACTOR*(b.x - a.x); // start with first guess
    c.f = fun(c.x);
    writeBracketToLog("findBracket init", bracket);

    while (b.f > c.f) { // loop until bracket found
        // compute parabolic extrapolation from a,b,c
        const double r = (b.x - a.x)*(b.f.value - c.f.value);
        const double q = (b.x - c.x)*(b.f.value - a.f.value);
        const double qmr = q - r; // used for sign
        const double qmr2 = 2.*std::max(fabs(qmr), MIN_DOUBLE); // prevent bad division
        h.x = b.x - ((b.x - c.x)*q - (b.x - a.x)*r)/(qmr >= 0. ? qmr2 : -qmr2);
        const double hlim = b.x + STRETCH_LIMIT*(c.x - b.x); // limit on hx

        // check the possible cases
        if ((b.x - h.x)*(h.x - c.x) > 0.) { // h is between b and c
            h.f = fun(h.x);
            if (h.f < c.f) { // minimum between b and c
                a = b;
                b = h;
                writeBracketToLog("findBracket final", bracket);
                return bracket;
            }
            else if (h.f > b.f) {// minimum between a and h
                c = h;
                writeBracketToLog("findBracket final", bracket);
                return bracket;
            }
            h.x = c.x + STRETCH_FACTOR*(c.x - b.x); // parabolic fit didn't work, use default stretch
            h.f = fun(h.x);
        }
        else if ((c.x - h.x)*(h.x - hlim) > 0.) { // fit is between c and hlim
            h.f = fun(h.x);
            if (h.f < c.f) {
                shift3(b.x, c.x, h.x, h.x + STRETCH_FACTOR*(h.x - c.x));
                shift3(b.f, c.f, h.f, fun(h.x));
            }
        }
        else if ((h.x - hlim)*(hlim - c.x) >= 0.) { // Limit h to max allowed value
            h.x = hlim;
            h.f = fun(h.x);
        }
        else { // reject h, use default stretch
            h.x = c.x + STRETCH_FACTOR*(c.x - b.x);
            h.f = fun(h.x);
        }
        // prepare next iteration and continue
        shift3(a.x, b.x, c.x, h.x);
        shift3(a.f, b.f, c.f, h.f);
        writeBracketToLog("findBracket step", bracket);
    }

    writeBracketToLog("findBracket final", bracket);
    return bracket;



/*    // bracket (+ shortcuts)
    NoisyBracket bracket {};
    NoisyIOPair1D & a = bracket.a;
    NoisyIOPair1D & b = bracket.b;
    NoisyIOPair1D & c = bracket.c;

    std::vector<double> newx(1); // helper array to invoke noisy function
    int count_newf = 0;

    a.x = initX;
    b.x = a.x + 0.25;

    if (++count_newf > MAX_NUM_EVAL_FOR_BRACKET) { _abortFindBracket(); }
    newx[0] = b.x;
    b.f = f1d.f(newx);
    writeBracketToLog("findBracket init", bracket);

    while (b.f == a.f) {  // if fb==fa, increase b until when the two values are different
        //cout << "WHILE 1    fa=" << a.f << "   fb=" << b.f << endl;
        b.x = a.x + 1.5*(b.x - a.x);
        //cout << "b=" << b.x << endl;

        if (++count_newf > MAX_NUM_EVAL_FOR_BRACKET) { _abortFindBracket(); }
        newx[0] = b.x;
        b.f = f1d.f(newx);
        writeBracketToLog("findBracket init", bracket);

        if (fabs(b.x) > hugeNum) {
            break;
        }
    }

    if (b.f == a.f) {  // could not find a b such that fb!=fa, now try looking in the other direction
        if (++count_newf > MAX_NUM_EVAL_FOR_BRACKET) { _abortFindBracket(); }

        b = a;
        a.x = b.x - 0.25;
        newx[0] = a.x;
        a.f = f1d.f(newx);
        writeBracketToLog("findBracket init", bracket);

        while (b.f == a.f) {   // if fb==fa, decrease a until when the two values are different
            //cout << "WHILE 1    fa=" << a.f << "   fb=" << b.f << endl;
            a.x = b.x - 1.5*(b.x - a.x);
            //cout << "a=" << a.x << endl;

            if (++count_newf > MAX_NUM_EVAL_FOR_BRACKET) { _abortFindBracket(); }
            newx[0] = a.x;
            a.f = f1d.f(newx);
            writeBracketToLog("findBracket init", bracket);

            if (fabs(a.x) > hugeNum) {
                throw std::runtime_error("[nfm::findBracket] Bracketing impossible. Cannot find a b such that fa!=fb");
            }
        }
    }

    if (b.f < a.f) {
        //cout << "COND b<a" << endl;
        c.x = b.x + (b.x - a.x);
        //cout << "c=" << c.x << endl;

        if (++count_newf > MAX_NUM_EVAL_FOR_BRACKET) { _abortFindBracket(); }
        newx[0] = c.x;
        c.f = f1d.f(newx);
        writeBracketToLog("findBracket (A)", bracket);

        while (!(c.f > b.f)) {
            //cout << "WHILE 2" << endl;
            if (c.f < b.f) {
                a = b;
                b = c;
                c.x = b.x + 2.*(b.x - a.x);
            }
            else {
                c.x = c.x + 2.*(c.x - b.x);
            }
            //cout << "c=" << c.x << endl;
            if (++count_newf > MAX_NUM_EVAL_FOR_BRACKET) { _abortFindBracket(); }
            newx[0] = c.x;
            c.f = f1d.f(newx);
            writeBracketToLog("findBracket (B)", bracket);
        }
    }
    else {
        c = b;
        b = a;
        //cout << "COND b>a" << endl;
        a.x = b.x - (c.x - b.x);
        //cout << "c=" << c.x << endl;

        if (++count_newf > MAX_NUM_EVAL_FOR_BRACKET) { _abortFindBracket(); }
        newx[0] = a.x;
        a.f = f1d.f(newx);
        writeBracketToLog("findBracket (C)", bracket);

        while (!(a.f > b.f)) {
            //cout << "WHILE 2" << endl;
            if (a.f < b.f) {
                c = b;
                b = a;
                a.x = b.x - 2.*(c.x - b.x);
            }
            else {
                a.x = a.x - 2.*(b.x - a.x);
            }
            //cout << "c=" << c.x << endl;

            if (++count_newf > MAX_NUM_EVAL_FOR_BRACKET) { abortFindBracket(); }
            newx[0] = a.x;
            a.f = f1d.f(newx);
            writeBracketToLog("findBracket (D)", bracket);
        }
    }

    writeBracketToLog("findBracket", bracket);
    */
}

NoisyIOPair1D brentMinimization(NoisyFunction &f1d, const double eps, NoisyBracket bracket)
{
    using namespace std;

    if (f1d.getNDim() != 1) {
        throw std::invalid_argument("[nfm::brentMinimization] The NoisyFunction is not 1D. Ndim=" + std::to_string(f1d.getNDim()));
    }



    /*
    int count = 0, lh = 0, rh = 0, cx = 0;
    bool force_lh = false, force_rh = false;

    if (b > a) { cout << "ERROR parabgoldMinimization(): b>a" << endl; }
    if (b > c) { cout << "ERROR parabgoldMinimization(): b>c" << endl; }

    double present_eps = (a > c)
                         ? (a.getF() - b.getF() - a.getDf() - b.getDf())
                         : (c.getF() - b.getF() - c.getDf() - b.getDf());

    NoisyValue x(1);
    double delta;
    double newf, dnewf;

    //cout << endl << endl << "a=" << a.getX(0) << "     b=" << b.getX(0) << "    c=" << c.getX(0) << endl;
    //cout << "fa=" << a.getF()  << "    fb=" << b.getF()  << "   fc=" << c.getF() << endl << endl;

    writeBracketToLog("parabgoldMinimization", a, b, c);

    while (present_eps > eps) {
        //update counter
        count++;
        //check the conditions for a,b and c
        if (b > a) { cout << "ERROR parabgoldMinimization(): b>a, count=" << count << endl; }
        if (b > c) { cout << "ERROR parabgoldMinimization(): b>c, count=" << count << endl; }

        //find the new x using a parabolic interpolation
        x.setX(b.getX(0) - 0.5*((b.getX(0) - a.getX(0))*(b.getX(0) - a.getX(0))*(b.getF() - c.getF()) - (b.getX(0) - c.getX(0))*(b.getX(0) - c.getX(0))*(b.getF() - a.getF()))/((b.getX(0) - a.getX(0))*(b.getF() - c.getF()) - (b.getX(0) - c.getX(0))*(b.getF() - a.getF())));

        //if the range has been reduced on the same side of b too many times,
        //or last time the parabolic prediction has led to a point equivalent to b,
        //use the golden bi-section rule on the other side of b
        if (((x.getX(0) < b.getX(0)) && (lh > 1)) || ((x.getX(0) > b.getX(0)) && (rh > 1)) || (cx > 0)) {
            if (b.getX(0) - a.getX(0) > c.getX(0) - b.getX(0)) {
                delta = (b.getX(0) - a.getX(0))*0.38197*0.5;
                x.setX(b.getX(0) - delta);
            }
            else {
                delta = (c.getX(0) - b.getX(0))*0.38197*0.5;
                x.setX(b.getX(0) + delta);
            }
        }

        //if it has been asked to propose a point on the left side of b
        if (force_lh) {
            delta = (b.getX(0) - a.getX(0))*0.38197*0.5;
            x.setX(b.getX(0) - delta);
            //cout << "force SX" << endl;
        }
        //if it has been asked to propose a point on the right side of b
        if (force_rh) {
            delta = (c.getX(0) - b.getX(0))*0.38197*0.5;
            x.setX(b.getX(0) + delta);
            //cout << "force DX" << endl;
        }

        //compute the value of f in the point x
        f1d.f(x.getX(), newf, dnewf);
        x.setF(newf, dnewf);

        //cout << "x=" << x.getX(0) << "     b=" << b.getX(0) << endl;
        force_rh = false;
        force_lh = false;
        if (x.getX(0) < b.getX(0)) {
            lh++;
            rh = 0;
            if (x > b) {
                cx = 0;
                //cout << "a=x" << endl;
                a = x;
            }
            else if (x < b) {
                cx = 0;
                //cout << "c=b" << endl;
                c = b;
                //cout << "b=x" << endl;
                b = x;
            }
            else {
                if (cx > 0) {
                    force_rh = true;
                    force_lh = false;
                }
                cx++;
                if (x.getF() < b.getF()) {
                    //cout << "b=x" << endl;
                    b = x;
                }
            }
        }
        else {
            lh = 0;
            rh++;
            if (x > b) {
                cx = 0;
                //cout << "c=x" << endl;
                c = x;
            }
            else if (x < b) {
                cx = 0;
                //cout << "a=b" << endl;
                a = b;
                //cout << "b=x" << endl;
                b = x;
            }
            else {
                if (cx > 0) {
                    force_lh = true;
                    force_rh = false;
                }
                cx++;
                if (x.getF() < b.getF()) {
                    //cout << "b=x" << endl;
                    b = x;
                }
            }
        }
        //update the actual eps value
        present_eps = (a > c)
                      ? (a.getF() - b.getF() - a.getDf() - b.getDf())
                      : (c.getF() - b.getF() - c.getDf() - b.getDf());
        //cout << "a=" << a.getX(0) << "     b=" << b.getX(0) << "    c=" << c.getX(0) << endl;
        //cout << "fa=" << a.getF()  << "    fb=" << b.getF()  << "   fc=" << c.getF() << endl << endl;
        if (cx > 2) { break; }
    }
    //cout << "ParabGold Terminated. count=" << count << endl;

    writeBracketToLog("parabgoldMinimization", a, b, c);
    */
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