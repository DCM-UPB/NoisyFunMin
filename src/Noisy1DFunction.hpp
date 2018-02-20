#ifndef NOISY_1D_FUN
#define NOISY_1D_FUN


class Noisy1DFunction
{
public:
    //Noisy 1D function
    virtual void f(const double &, double &, double &) = 0;
    //              ^input          ^output   ^error on the output
};


#endif
