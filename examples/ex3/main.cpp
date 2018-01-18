#include <iostream>
#include <cmath>
#include <math.h>
#include <fstream>
#include <random>

#include "NoisyFunction.hpp"
#include "ConjGrad.hpp"

#include "FeedForwardNeuralNetwork.hpp"
#include "PrintUtilities.hpp"


// exp(-a*(x-b)^2)
class Gaussian{
private:
   double _a;
   double _b;
   
public:
   Gaussian(a, b){
      _a = a;
      _b = b;
   }
   
   double f(const double &x){
      return exp(-a*pow(x-b, 2));
   }
};



class FitNN: public NoisyFunctionWithGradient{
private:
   FeedForwardNeuralNetwork * _ffnn;
   Gaussian * _gss;
   
public:
   FitNN(FeedForwardNeuralNetwork * ffnn, Gaussian * gss): NoisyFunctionWithGradient(ffnn->getNBeta()){
      _ffnn = ffnn;
      _gss = gss;
   }
   
   void f(const double * in, double &f, double &df){
      const double x = 0.;
      
      _ffnn->setInput(&x, 1);
      for (int i=0; i<_ffnn->getNBeta(); ++i){
         _ffnn->setBeta(i, in[i]);
      }
      _ffnn->propagateFFNN();
      f = pow(_ffnn->getOutput(1) - _gss->(x), 2);
      df = 0.;
   }
   
   void grad(const double * in, double * g, double * dg){
      for (int i=0; i<_ffnn->getNBeta(); ++i){
         _ffnn->setBeta(i, in[i]);
      }
      _ffnn->propagateFFNN();
      for (int i=0; i<_ffnn->getNBeta(); ++i){
         g[i] = _ffnn->getVariationalFirstDerivative(1, i);
         dg[i] = 0.;
      }
   }
};





int main() {
    using namespace std;
    
    FeedForwardNeuralNetwork * ffnn = new FeedForwardNeuralNetwork(...);
    Gaussian * gauss = new Gaussian(...);
    
    FitNN * fitNN = new FitNN(ffnn, gauss);
    
    ... conj->setX(...);  // read from the ffnn:    ffnn->getBeta()
    
    ... conj->setEpsX(0.01);
    
    ConjGrad * conj = new ConjGrad(fitNN);
    conj->findMin();
    
    ... conj->getX(...); // you will get the "right" beta
    ... ffnn->setBeta(...);
    
    delete conj;
    delete fitNN;
    delete gauss;
    delete ffnn;

    // end
    return 0;
}
