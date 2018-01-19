#include <iostream>
#include <cmath>
#include <math.h>
#include <fstream>
#include <random>

#include "NoisyFunction.hpp"
#include "ConjGrad.hpp"
#include "MCIntegrator.hpp"
#include "FeedForwardNeuralNetwork.hpp"
#include "PrintUtilities.hpp"


//1D Target Function
class Function1D {
public:
  //Function1D(){}
  virtual ~Function1D(){}

  // Function evaluation
  virtual double f(const double &) = 0;
  //                    ^input
};

// exp(-a*(x-b)^2)
class Gaussian: public Function1D{
private:
  double _a;
  double _b;

public:
  Gaussian(const double &a, const double &b): Function1D(){
    _a = a;
    _b = b;
  }

  double f(const double &x){
    return exp(-_a*pow(x-_b, 2));
  }
};

// Quadratic distance between NN and 1D target function
class NNFitDistance1D: public MCIObservableFunctionInterface{
private:
  FeedForwardNeuralNetwork * _ffnn;
  Function1D * _ftarget;
public:
  NNFitDistance1D(FeedForwardNeuralNetwork * ffnn, Function1D * ftarget): MCIObservableFunctionInterface(1, 1) {_ffnn = ffnn; _ftarget = ftarget;}
protected:
    void observableFunction(const double * in, double * out){

      _ffnn->setInput(1, &in[0]);
      _ffnn->FFPropagate();
      out[0] = pow(_ffnn->getOutput(1) - _ftarget->f(in[0]), 2);

    }
};

// Gradient of Quadratic distance between NN and 1D target function
class NNFitGradient1D: public MCIObservableFunctionInterface{
private:
  FeedForwardNeuralNetwork * _ffnn;
  Function1D * _ftarget;
public:
  NNFitGradient1D(FeedForwardNeuralNetwork * ffnn, Function1D * ftarget): MCIObservableFunctionInterface(1, ffnn->getNBeta()) {_ffnn = ffnn; _ftarget = ftarget;}
protected:
  void observableFunction(const double * in, double * out){

    _ffnn->setInput(1, &in[0]);
    _ffnn->FFPropagate();

    double nnout = _ffnn->getOutput(1);
    double fx = _ftarget->f(in[0]);

    for(int i=0; i<_nobs; ++i){
        out[i] = 2.*(nnout - fx)*_ffnn->getVariationalFirstDerivative(1, i);
      }
  }
};

//
class FitNN1D: public NoisyFunctionWithGradient{
private:
  FeedForwardNeuralNetwork * _ffnn;
  Function1D * _ftarget;
  MCI * _mcif, * _mcig;
  long _nmc;
  double _x0, _step0;

public:
  FitNN1D(FeedForwardNeuralNetwork * ffnn, Function1D * ftarget, const long &nmc, double * irange): NoisyFunctionWithGradient(ffnn->getNBeta()){
    _ffnn = ffnn;
    _ftarget = ftarget;
    _nmc = nmc;

    _x0 = 0.5*(irange[0] + irange[1]);
    _step0 = 0.1*(irange[1] - irange[0]);

    _mcif = new MCI(1);
    _mcig = new MCI(1);

    MCIObservableFunctionInterface * fobs = new NNFitDistance1D(_ffnn, _ftarget);
    MCIObservableFunctionInterface * gobs = new NNFitGradient1D(_ffnn, _ftarget);

    _mcif->addObservable(fobs);
    _mcig->addObservable(gobs);

    double ** intrange = new double*[1];
    intrange[0] = irange;
    _mcif->setIRange(intrange);
    _mcig->setIRange(intrange);

    double targetacceptrate = 0.5;
    _mcif->setTargetAcceptanceRate(&targetacceptrate);
    _mcig->setTargetAcceptanceRate(&targetacceptrate);

  }

  void f(const double * in, double &f, double &df){
    int i;

    //set new NN betas
    for (i=0; i<_ffnn->getNBeta(); ++i){
      _ffnn->setBeta(i, in[i]);
    }

    // initial walker position and step size
    _mcif->setX(&_x0);
    _mcif->setMRT2Step(&_step0);

    // integrate
    double * average = new double[_mcif->getNObsDim()];
    double * error = new double[_mcif->getNObsDim()];
    _mcif->integrate(_nmc, average, error);

    // write out results
    f = average[0];
    df = error[0];
  }

  void grad(const double * in, double * g, double * dg){
    int i, nbeta = _ffnn->getNBeta();

    //set new NN betas
    for (i=0; i<nbeta; ++i){
      _ffnn->setBeta(i, in[i]);
    }

    // initial walker position and step size
    _mcig->setX(&_x0);
    _mcig->setMRT2Step(&_step0);

    // integrate
    double * average = new double[_mcig->getNObsDim()];
    double * error = new double[_mcig->getNObsDim()];
    _mcig->integrate(_nmc, average, error);

    // write out results
    for(i=0; i<nbeta; ++i){
      g[i] = average[i];
      dg[i] = error[i];
    }
  }

};


int main() {
  using namespace std;

  int nl, nh1, nh2;
  long nmc = 40000l;
  //double epsx = 0.01;
  double * irange = new double[2];
  irange[0] = -10.;
  irange[1] = 10.;

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
  cout << "This means creating the target function, the quadratic distance NoisyFunction to optimize and the Conjugate Gradient object." << endl << endl;
  cin.ignore();

  // NON I/O CODE
  Gaussian * gauss = new Gaussian(0.1,0.);
  FitNN1D * fitnn = new FitNN1D(ffnn, gauss, nmc, irange);
  ConjGrad * conj = new ConjGrad(fitnn);
  //

  cout << "First copy the initial betas from the NN to CG and set epsx." << endl << endl;
  cin.ignore();

  // NON I/O CODE
  for(int i = 0; i<ffnn->getNBeta(); ++i) conj->setX(i, ffnn->getBeta(i));
  //  conj->setEpsX(epsx);
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

  cout << "Let us compare the NN to target function values:" << endl;
  cin.ignore();

  double x=-10., dx=1.;
  for(int i=0; i<21; ++i) {
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
