#include "nfm/1DTools.hpp"

#include "nfm/LogManager.hpp"

#include <string>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <iostream>

const double hugeNum = 1000000.;
const int MAX_NUM_EVAL_FOR_BRACKET = 100;

void _abortFindBracket()
{
    throw std::runtime_error("NoisyFunctionMin Error! Bracketing is taking way too long...");
}

namespace nfm
{

void writeBracketToLog(const std::string &key, const NoisyBracket &bracket)
{
    using namespace std;

    stringstream s;
    s << key << ":    " <<
      bracket.a.x << " -> " << bracket.a.f << "    " <<
      bracket.b.x << " -> " << bracket.b.f << "    " <<
      bracket.c.x << " -> " << bracket.c.f;
    s << flush;
    LogManager::logString(s.str(), LogLevel::VERBOSE);
}

NoisyBracket findBracket(NoisyFunction &f1d, const double initX)
{
    using namespace std;

    if (f1d.getNDim() != 1) {
        throw std::invalid_argument("[nfm::findBracket] The NoisyFunction is not 1D. Ndim=" + to_string(f1d.getNDim()));
    }

    // bracket (+ shortcuts)
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

            if (++count_newf > MAX_NUM_EVAL_FOR_BRACKET) { _abortFindBracket(); }
            newx[0] = a.x;
            a.f = f1d.f(newx);
            writeBracketToLog("findBracket (D)", bracket);
        }
    }

    writeBracketToLog("findBracket", bracket);
}


void parabgoldMinimization(NoisyFunction &f1d, const double eps, NoisyValue &a, NoisyValue &b, NoisyValue &c)
{
    using namespace std;

    if (f1d.getNDim() != 1) {
        throw std::invalid_argument("[nfm::parabgoldMinimization] The NoisyFunction is not 1D. Ndim=" + std::to_string(f1d.getNDim()));
    }

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
}
} // namespace nfm