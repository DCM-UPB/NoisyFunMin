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

    auto * log_manager = new LogManager();
    log_manager->setLoggingOn();
    log_manager->setLoggingPathFile("log.txt");

    double x1 = 0., x2 = 3.;
    double f1 = 5., f2 = 4., df1 = 1.1, df2 = 0.2;
    NoisyValue p1(1);
    p1.setX(x1);
    p1.setF(f1, df1);
    NoisyValue p2(1);
    p2.setX(x2);
    p2.setF(f2, df2);

    cout << " - - - Check NoisyValue" << endl;
    cout << "f1<f2 ? (0 expected) " << (p1 < p2) << endl;
    cout << "f1<=f2 ? (1 expected) " << (p1 <= p2) << endl;
    cout << "f1>f2 ? (0 expected) " << (p1 > p2) << endl;
    cout << "f1>=f2 ? (1 expected) " << (p1 >= p2) << endl;
    cout << "f1==f2 ? (1 expected) " << (p1 == p2) << endl << endl;

    //check bracketing
    cout << " - - - Check nfm::findBracket()" << endl;
    NoisyValue p3(1);
    p1.setX(10.1);
    auto * f1d = new F1D();
    f1d->f(p1.getX(), f1, df1);
    p1.setF(f1, df1);
    nfm::findBracket(f1d, p1, p2, p3);
    cout << "a=" << p1.getX(0) << "     b=" << p2.getX(0) << "     c=" << p3.getX(0) << endl;
    cout << "fa=" << p1.getF() << "      fb=" << p2.getF() << "      fc=" << p3.getF() << endl << endl;
    log_manager->logString("\n\n=========================================================================\n\n");

    // check parabgold
    cout << " - - - Check nfm::parabgoldMinimization()" << endl;
    nfm::parabgoldMinimization(f1d, 0., p1, p2, p3);
    cout << "Minimum of f1d is " << p2.getF() << " +- " << p2.getDf() << "    in " << p2.getX(0) << endl << endl;
    log_manager->logString("\n\n=========================================================================\n\n");

    // check Conjugate Gradient
    cout << " - - - Check ConjGrad" << endl;
    F3D * f3d = new F3D();
    ConjGrad cjgrad(f3d);
    double x[3];
    x[0] = -2.;
    x[1] = 1.0;
    x[2] = 0.0;
    cjgrad.setX(x);
    cjgrad.findMin();
    cout << "Minimum of f3d is in " << cjgrad.getX(0) << "   " << cjgrad.getX(1) << "   " << cjgrad.getX(2) << endl;
    cout << "Value of the minimum is " << cjgrad.getF() << " +- " << cjgrad.getDf() << endl << endl;

    delete f3d;
    delete f1d;
    delete log_manager;

    return 0;
}
