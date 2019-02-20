#ifndef TEST_NFM_FUNCTIONS
#define TEST_NFM_FUNCTIONS

#include <cmath>
#include "nfm/NoisyFunction.hpp"


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


#endif
