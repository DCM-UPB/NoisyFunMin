#include <iostream>
#include <assert.h>

#include "NoisyFunction.hpp"
#include "NoisyFunctionValue.hpp"
#include "1DTools.hpp"
#include "ConjGrad.hpp"
#include "LogNFM.hpp"




class Parabola: public NoisyFunction
{
public:
    Parabola():NoisyFunction(1){}

    void f(const double * in, double &f, double &df)
    {
        f=in[0]*in[0];
        df=0.5;
    }
};



class Well: public NoisyFunction
{
public:
    Well():NoisyFunction(1){}

    void f(const double * in, double &f, double &df)
    {
        if ((in[0] <= -1.) || (in[0] >= 1.)){
            f = 1.;
        } else {
            f = -1.;
        }
        df=0.1;
    }
};




int main(){

    using namespace std;

    NFMLogManager * log_manager = new NFMLogManager();

    double f, df;

    // define 3 noisy function values
    NoisyFunctionValue p1(1);
    NoisyFunctionValue p2(1);
    NoisyFunctionValue p3(1);



    // check parabola   x^2   ...
    Parabola parabola;

    // ... starting from x=-1000
    p1.setX(-1000.);
    parabola.f(p1.getX(), f, df);
    p1.setF(f, df);
    nfm::findBracket(&parabola, p1, p2, p3);
    assert(p1.getX(0) < 0.);
    assert(p3.getX(0) > 0.);

    // ... starting from x=1000
    log_manager->writeOnLog("\n\n=========================================================================\n\n");
    p1.setX(1000.);
    parabola.f(p1.getX(), f, df);
    p1.setF(f, df);
    nfm::findBracket(&parabola, p1, p2, p3);
    assert(p1.getX(0) < 0.);
    assert(p3.getX(0) > 0.);

    // ... starting from x=0
    log_manager->writeOnLog("\n\n=========================================================================\n\n");
    p1.setX(0.);
    parabola.f(p1.getX(), f, df);
    p1.setF(f, df);
    nfm::findBracket(&parabola, p1, p2, p3);
    assert(p1.getX(0) < 0.);
    assert(p3.getX(0) > 0.);



    // check well function   -1 if (-1 < x < 1) else +1   ...
    Well well;

    // ... starting from x=-10
    log_manager->writeOnLog("\n\n=========================================================================\n\n");
    p1.setX(-10.);
    well.f(p1.getX(), f, df);
    p1.setF(f, df);
    nfm::findBracket(&well, p1, p2, p3);
    assert(p1.getX(0) < -1.);
    assert(p3.getX(0) > 1.);

    // ... starting from x=-1000
    log_manager->writeOnLog("\n\n=========================================================================\n\n");
    p1.setX(-1000.);
    well.f(p1.getX(), f, df);
    p1.setF(f, df);
    bool flag_exception_thrown = false;
    try{
        nfm::findBracket(&well, p1, p2, p3);
    } catch(exception& e){
        flag_exception_thrown = true;
    }
    assert(flag_exception_thrown);

    // ... starting from x=10
    log_manager->writeOnLog("\n\n=========================================================================\n\n");
    p1.setX(10.);
    well.f(p1.getX(), f, df);
    p1.setF(f, df);
    nfm::findBracket(&well, p1, p2, p3);
    assert(p1.getX(0) < -1.);
    assert(p3.getX(0) > 1.);


    return 0;
}
