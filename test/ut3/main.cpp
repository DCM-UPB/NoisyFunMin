#include <iostream>
#include <assert.h>
#include <math.h>

#include "nfm/NoisyFunction.hpp"
#include "nfm/NoisyFunctionValue.hpp"
#include "nfm/1DTools.hpp"
#include "nfm/ConjGrad.hpp"
#include "nfm/DynamicDescent.hpp"
#include "nfm/LogNFM.hpp"



class F3D: public NoisyFunctionWithGradient
{
public:
    F3D():NoisyFunctionWithGradient(3){}

    void f(const double * in, double &f, double &df)   // f = (x-1)^4 + (y+1.5)^4 + (z-0.5)^4
    {
        f=pow(in[0]-1.,4)+pow(in[1]+1.5,4)+pow(in[2]-0.5,4);
        df=0.00001;
    }

    void grad(const double *in, double *g, double *dg)
    {
        g[0]=4.*pow( in[0]-1.0, 3);
        g[1]=4.*pow( in[1]+1.5, 3);
        g[2]=4.*pow( in[2]-0.5, 3);
        dg[0]=0.000001; dg[1]=0.000001; dg[2]=0.000001;
    }

};


int main(){
    using namespace std;

    NFMLogManager * log_manager = new NFMLogManager();

    // define 3D function that I want to minimise
    F3D * f3d = new F3D();
    // introduce array with the initial position
    double x[3];


    // test ConjGrad
    ConjGrad cjgrad(f3d);
    x[0] = -2.;   x[1] = 1.0;   x[2] = 0.0;
    cjgrad.setX(x);
    cjgrad.findMin();
    assert(std::abs(cjgrad.getX(0)-1.0) < 0.1);
    assert(std::abs(cjgrad.getX(1)+1.5) < 0.1);
    assert(std::abs(cjgrad.getX(2)-0.5) < 0.1);

    delete f3d;
    delete log_manager;

    return 0;
}
