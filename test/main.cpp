#include <cmath>
#include <iostream>

#include "nfm/1DTools.hpp"
#include "nfm/ConjGrad.hpp"
#include "nfm/LogManager.hpp"

#include "TestNFMFunctions.hpp"


int main()
{
    using namespace std;
    using namespace nfm;

    LogManager::setLogLevel(LogLevel::VERBOSE);
    //LogManager::setLoggingPathFile("log.txt");

    double x1 = 0., x2 = 3.;
    double f1 = 5., f2 = 4., df1 = 1.1, df2 = 0.2;
    NoisyIOPair1D p1{x1, NoisyValue({f1, df1})};
    NoisyIOPair1D p2{x2, NoisyValue({f2, df2})};
    vector<double> xvec(1); // 1d vector to call NoisyFunctions

    cout << " - - - Check NoisyValue" << endl;
    cout << "f1<f2 ? (0 expected) " << (p1.f < p2.f) << endl;
    cout << "f1<=f2 ? (1 expected) " << (p1.f <= p2.f) << endl;
    cout << "f1>f2 ? (0 expected) " << (p1.f > p2.f) << endl;
    cout << "f1>=f2 ? (1 expected) " << (p1.f >= p2.f) << endl;
    cout << "f1==f2 ? (1 expected) " << (p1.f == p2.f) << endl << endl;

    //check bracketing
    cout << " - - - Check nfm::findBracket()" << endl;
    F1D f1d;
    p1.x = 10.1;
    xvec[0] = p1.x;
    p1.f = f1d.f(xvec);
    NoisyBracket brk = nfm::findBracket(f1d, p1.x, p2.x);
    cout << "a=" << brk.a.x << "     b=" << brk.b.x << "     c=" << brk.c.x << endl;
    cout << "fa=" << brk.a.f << "      fb=" << brk.b.f << "      fc=" << brk.c.f << endl << endl;
    LogManager::logString("\n\n=========================================================================\n\n");

    // check parabgold
    cout << " - - - Check nfm::brentMinimization()" << endl;
    NoisyIOPair1D p3 = nfm::brentMinimization(f1d, brk, 0.);
    cout << "Minimum of f1d is " << p3.f << "    in " << p3.x << endl << endl;
    LogManager::logString("\n\n=========================================================================\n\n");

    // check Conjugate Gradient
    cout << " - - - Check ConjGrad" << endl;
    F3D f3d;
    ConjGrad cjgrad(&f3d);
    double x[3];
    x[0] = -2.;
    x[1] = 1.0;
    x[2] = 0.0;
    cjgrad.setX(x);
    cjgrad.setEpsF(0.01);
    cjgrad.findMin();
    cout << "Minimum of f3d is in " << cjgrad.getX(0) << "   " << cjgrad.getX(1) << "   " << cjgrad.getX(2) << endl;
    cout << "Value of the minimum is " << cjgrad.getF() << " +- " << cjgrad.getDf() << endl << endl;

    return 0;
}
