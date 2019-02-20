#include <iostream>
#include <assert.h>

#include "nfm/NoisyFunction.hpp"
#include "nfm/NoisyFunctionValue.hpp"
#include "nfm/1DTools.hpp"
#include "nfm/ConjGrad.hpp"
#include "nfm/LogNFM.hpp"




class PowerFour: public NoisyFunction
{
public:
    PowerFour():NoisyFunction(1){}

    void f(const double * in, double &f, double &df)
    {
        f=in[0]*in[0]*in[0]*in[0];
        df=0.000000001;
    }
};






int main(){
    using namespace std;

    NFMLogManager log_manager;
    //log_manager.setLoggingOn();

    double f, df;

    // define 3 noisy function values
    NoisyFunctionValue p1(1);
    NoisyFunctionValue p2(1);
    NoisyFunctionValue p3(1);



    // check pwr4   x^4   ...
    PowerFour pwr4;

    // ... using a=-3.   b=-2.   c=5.
    p1.setX(-3.);   pwr4.f(p1.getX(), f, df);   p1.setF(f, df);
    p2.setX(-2.);   pwr4.f(p2.getX(), f, df);   p2.setF(f, df);
    p3.setX(+5.);   pwr4.f(p3.getX(), f, df);   p3.setF(f, df);
    nfm::parabgoldMinimization(&pwr4, 0., p1, p2, p3);
    assert(p2.getX(0)<0.1);
    assert(p2.getX(0)>-0.1);
    assert(p2.getF()<0.00001);


    return 0;
}
