#include <iostream>
#include <cmath>

#include "NoisyFunction.hpp"
#include "NoisyFunctionValue.hpp"
#include "1DTools.hpp"
#include "ConjGrad.hpp"
#include "LogNFM.hpp"


class F1D: public NoisyFunction
{
public:
    F1D():NoisyFunction(1){}

    void f(const double * in, double &f, double &df)
    {
        f=pow(*in-1.,4);
        df=0.00001;
    }
};


class F3D: public NoisyFunctionWithGradient
{
public:
    F3D():NoisyFunctionWithGradient(3){}

    void f(const double * in, double &f, double &df)
    {
        f=pow(in[0]-1.,4)+pow(in[1]+1.5,4)+pow(in[2]-0.5,4);
        df=0.00001;
    }

    void grad(const double *in, double *g, double *dg)
    {
        g[0]=4.*pow(in[0]-1.,3);
        g[1]=4.*pow(in[1]+1.5,3);
        g[2]=4.*pow(in[2]-0.5,3);
        dg[0]=0.000001; dg[1]=0.000001; dg[2]=0.000001;
    }

};


int main(){
    using namespace std;

    NFMLogManager * log_manager = new NFMLogManager();
    log_manager->setLoggingOn();
    log_manager->setLoggingPathFile("log.txt");

    double x1=0. , x2=3.;
    double f1=5. , f2=4. , df1=1.1 , df2=0.2;
    NoisyFunctionValue p1(1);
    p1.setX(x1);
    p1.setF(f1,df1);
    NoisyFunctionValue p2(1);
    p2.setX(x2);
    p2.setF(f2,df2);

    cout << " - - - Check NoisyFunctionValue" << endl;
    cout << "f1<f2 ? (0 expected) " << (p1<p2) << endl;
    cout << "f1<=f2 ? (1 expected) " << (p1<=p2) << endl;
    cout << "f1>f2 ? (0 expected) " << (p1>p2) << endl;
    cout << "f1>=f2 ? (1 expected) " << (p1>=p2) << endl;
    cout << "f1==f2 ? (1 expected) " << (p1==p2) << endl << endl;

    //check bracketing
    cout << " - - - Check nfm::findBracket()" << endl;
    NoisyFunctionValue p3(1);
    p1.setX(10.1);
    F1D * f1d = new F1D();
    f1d->f(p1.getX(), f1, df1);
    p1.setF(f1,df1);
    nfm::findBracket(f1d, p1, p2, p3);
    cout << "a="  << p1.getX(0) << "     b="  << p2.getX(0) << "     c=" << p3.getX(0) << endl;
    cout << "fa=" << p1.getF() << "      fb=" << p2.getF() << "      fc=" << p3.getF() << endl << endl;
    log_manager->writeOnLog("\n\n=========================================================================\n\n");

    // check parabgold
    cout << " - - - Check nfm::parabgoldMinimization()" << endl;
    nfm::parabgoldMinimization(f1d, 0., p1, p2, p3);
    cout << "Minimum of f1d is " << p2.getF() << " +- " << p2.getDf() << "    in " << p2.getX(0) << endl << endl;
    log_manager->writeOnLog("\n\n=========================================================================\n\n");

    // check Conjugate Gradient
    cout << " - - - Check ConjGrad" << endl;
    F3D * f3d = new F3D();
    ConjGrad cjgrad(f3d);
    double x[3];
    x[0] = -2.;   x[1] = 1.0;   x[2] = 0.0;
    cjgrad.setX(x);
    cjgrad.findMin();
    cout << "Minimum of f3d is in " << cjgrad.getX(0) << "   " << cjgrad.getX(1) << "   " << cjgrad.getX(2) << endl;
    cout << "Value of the minimum is " << cjgrad.getF() << " +- " << cjgrad.getDf() << endl << endl;

    delete f3d;
    delete f1d;
    delete log_manager;

    return 0;
}
