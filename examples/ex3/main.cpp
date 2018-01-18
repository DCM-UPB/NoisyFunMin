#include <iostream>
#include <cmath>
#include <math.h>
#include <fstream>
#include <random>

#include "NoisyFunction.hpp"
#include "ConjGrad.hpp"

#include "FeedForwardNeuralNetwork.hpp"
#include "PrintUtilities.hpp"

//1D Target Function
class Function1D {
public:
  //Function1D(){}
  virtual ~Function1D(){}

  // Function evaluation
  virtual double f(const double) = 0;
  //                    ^input
};


// exp(-a*(x-b)^2)
class Gaussian: public Function1D{
private:
  double _a;
  double _b;

public:
  Gaussian(const double a, const double b): Function1D(){
    _a = a;
    _b = b;
  }

  double f(const double x){
    return exp(-_a*pow(x-_b, 2));
  }
};


//Takes
class FitNN1D: public NoisyFunctionWithGradient{
private:
  FeedForwardNeuralNetwork * _ffnn;
  Function1D * _ftarget;
  double _nmc;

public:
  FitNN1D(FeedForwardNeuralNetwork * ffnn, Function1D * ftarget, const int nmc): NoisyFunctionWithGradient(ffnn->getNBeta()){
    _ffnn = ffnn;
    _ftarget = ftarget;
    _nmc = nmc;
  }

  void f(const double * in, double &f, double &df){
    int i;
    double x=-10., dx=20./_nmc;

    for (i=0; i<_ffnn->getNBeta(); ++i){
      _ffnn->setBeta(i, in[i]);
    }

    f = 0.;
    df = 0.;
    for(i=0; i<_nmc; ++i){
      _ffnn->setInput(1, &x);
      x+=dx;
      _ffnn->FFPropagate();
      f += pow(_ffnn->getOutput(1) - _ftarget->f(x), 2);
    }
    f /= _nmc;
  }

  void grad(const double * in, double * g, double * dg){
    int i, j, nbeta = _ffnn->getNBeta();
    double x=-10., dx=20./_nmc, nnout, fx;

    for (j=0; j<nbeta; ++j){
      _ffnn->setBeta(j, in[j]);
      g[j] = 0.;
      dg[j] = 0.;
    }

    for(i=0; i<_nmc; ++i){
      _ffnn->setInput(1, &x);
      x+=dx;
      _ffnn->FFPropagate();
      nnout = _ffnn->getOutput(1);
      fx = _ftarget->f(x);
      for(j=0; j<nbeta; ++j) {
        g[j] += 2.*(nnout - fx)*_ffnn->getVariationalFirstDerivative(1, j);
      }
    }
    for (j=0; j<nbeta; ++j){
      g[j] /= _nmc;
    }
  }

};


int main() {
  using namespace std;

  int nl, nh1, nh2;
  int nmc = 200;
  double epsx = 0.01;

  cout << "Let's start by creating a Feed Forward Artificial Neural Netowrk (FFANN)" << endl;
  cout << "========================================================================" << endl;
  cin.ignore();

  cout << "How many units should the first hidden layer(s) have? ";
  cin >> nh1;
  cout << "How many units should the second hidden layer(s) have? (<=1 for no second hidden layer)";
  cin >> nh2;

  nl = (nh2>1)? 4:3; 
  cout << "We generate a FFANN with " << nl << " layers and 2, " << nh1;
  if (nh2>1) { cout << ", " << nh2;}
  cout << ", 2 units respectively" << endl;
  cout << "========================================================================" << endl << endl;
  cin.ignore();


  // NON I/O CODE
  FeedForwardNeuralNetwork * ffnn = new FeedForwardNeuralNetwork(2, nh1, 2);
  if (nh2>1) ffnn->pushHiddenLayer(nh2);
  //


  cout << "Graphically it looks like this" << endl;
  cin.ignore();
  printFFNNStructure(ffnn);
  cout << endl << endl;
  cin.ignore();

  cout << "Connecting the FFNN..." << endl;
  cout << endl << endl;
  cin.ignore();

  // NON I/O CODE
  ffnn->connectFFNN();
  //

  cout << "Adding derivatives substrates..." << endl;
  cout << endl << endl;
  cin.ignore();

  // NON I/O CODE
  //ffnn->addFirstDerivativeSubstrate();
  //ffnn->addSecondDerivativeSubstrate();
  ffnn->addVariationalFirstDerivativeSubstrate();
  //

  cout << "The Neural Network is now prepared." << endl << endl;
  cin.ignore();

  cout << "Now we set up the Noisy Function Minimization via Conjugate Gradient." << endl;
  cout << "This means creating the target function, the quadratic distance NoisyFunction to optimize and the Conjuagte Gradient object." << endl << endl;
  cin.ignore();

  // NON I/O CODE
  Gaussian * gauss = new Gaussian(0.01,0.);
  FitNN1D * fitnn = new FitNN1D(ffnn, gauss, nmc);
  ConjGrad * conj = new ConjGrad(fitnn);
  //

  cout << "First copy the initial betas from the NN to CG and set epsx." << endl << endl;
  cin.ignore();

  // NON I/O CODE
  for(int i = 0; i<ffnn->getNBeta(); ++i) conj->setX(i, ffnn->getBeta(i));
  conj->setEpsX(epsx);
  //

  cout << "Let CG find the optimal betas..." << endl << endl;
  cin.ignore();

  // NON I/O CODE
  conj->findMin();
  //

  cout << "Done." << endl;
  cout << "Now copy over the optimal betas from CG to the NN." << endl << endl;
  cin.ignore();

  for(int i = 0; i<ffnn->getNBeta(); ++i) ffnn->setBeta(i, conj->getX(i));

  cout << "Finally the NN is fitted!" << endl << endl;
  cin.ignore();

  cout << "Let's print the betas:" << endl;
  cin.ignore();
  for(int i = 0; i<ffnn->getNBeta(); ++i) cout << i << ": " << ffnn->getBeta(i) << endl;;

  cout << endl;

  cout << "And now we compare NN to target function values:" << endl;
  cin.ignore();

  double x=-10., dx=20./nmc;
  for(int i=0; i<nmc; ++i) {
    ffnn->setInput(1, &x);
    ffnn->FFPropagate();
    cout << "x: " << x << " f(x): " << gauss->f(x) << " nn(x): " << ffnn->getOutput(1) << endl;
    x+=dx;
  }
  cout << endl;


  delete conj;
  delete fitnn;
  delete gauss;
  delete ffnn;

  // end
  return 0;
}
