#include <iostream>
#include <cmath>
#include <math.h>
#include <fstream>
#include <random>

#include "nfm/NoisyFunction.hpp"
#include "nfm/Adam.hpp"
#include "nfm/LogNFM.hpp"



class Noiseless2DParabola: public NoisyFunctionWithGradient{
public:
    Noiseless2DParabola(): NoisyFunctionWithGradient(2){}

    void f(const double * in, double &f, double &df){
        f = pow(in[0]-1., 2) + pow(in[1]+2., 2);   // minimum in (1, -2)
        df = 0.0;
    }

    void grad(const double * in, double * g, double * dg){
        g[0] = 2. * (in[0] - 1.);
        g[1] = 2. * (in[1] + 2.);
        dg[0] = 0.0;
        dg[1] = 0.0;
    }
};



class Noisy2DParabola: public NoisyFunctionWithGradient{
private:
    const double _sigma = 0.15;
    std::random_device _rdev;
    std::mt19937_64 _rgen;
    std::uniform_real_distribution<double> _rd;  //after initialization (done in the constructor) can be used with _rd(_rgen)

public:
    Noisy2DParabola(): NoisyFunctionWithGradient(2){
        // initialize random generator
        _rgen = std::mt19937_64(_rdev());
        _rd = std::uniform_real_distribution<double>(-_sigma, _sigma);
    }

    void f(const double * in, double &f, double &df){
        f = pow(in[0]-1., 2) + pow(in[1]+2., 2);   // minimum in (1, -2)
        df = _sigma;
        f += _rd(_rgen);
    }

    void grad(const double * in, double * g, double * dg){
        g[0] = 2. * (in[0] - 1.);
        g[1] = 2. * (in[1] + 2.);
        dg[0] = 2.*_sigma;
        dg[1] = 2.*_sigma;
        g[0] += 2.*_rd(_rgen);
        g[1] += 2.*_rd(_rgen);
    }
};




int main() {
    using namespace std;

    cout << "We want to minimize the 2D function" << endl;
    cout << "    (x-1)^2 + (y+2)^2" << endl;
    cout << "whose min is in (1, -2)." << endl << endl << endl;

    NFMLogManager log;        
    //log.setLoggingOn(); // use this to enable log printout
    //log.setLogLevel(2); // use this for verbose printout of the Adam method

    cout << "we first minimize it, supposing to have no noise at all" << endl;

    Noiseless2DParabola * nlp = new Noiseless2DParabola();

    Adam * adam = new Adam(nlp);

    double initpos[2];
    initpos[0] = -1.;
    initpos[1] = -1.;
    adam->setX(initpos);

    adam->findMin();

    cout << "The found minimum is: ";
    cout << adam->getX(0) << "    " << adam->getX(1) << endl << endl << endl;



    cout << "Now we repeat the minimisation adding a noise to the function and its gradient." << endl;

    delete adam;
    delete nlp;

    Noisy2DParabola * np = new Noisy2DParabola();

    // in noisy-target low-dim cases like this we often need to change adam parameters from default (in brackets)
    adam = new Adam(np, true /* calculate/use gradient error for stopping */, 20 /* max n constant (within error) values before stopping */,
                    true /* use averaging to calculate final parameters */, 1.0 /* step size factor (0.001) */);
    adam->setX(initpos);

    adam->findMin();

    cout << "The found minimum is: ";
    cout << adam->getX(0) << "    " << adam->getX(1) << endl << endl;


    delete np;
    delete adam;


    // end
    return 0;
}
