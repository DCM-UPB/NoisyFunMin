#include <iostream>
#include <cmath>
#include <math.h>
#include <pthread.h>

#include "NoisyFunction.hpp"
#include "ConjGrad.hpp"
#include "MCIntegrator.hpp"
#include "FeedForwardNeuralNetwork.hpp"
#include "PrintUtilities.hpp"
// thread0.c
#include <pthread.h>



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

// Sampling function 1D depending on Function1D object, here uniform placeholder
// It's not necessary to define unifrom distributions explicitly, but done here for demonstration purposes
class SamplingFunction1D: public MCISamplingFunctionInterface{
private:
  Function1D * _f1d;
public:
  SamplingFunction1D(Function1D * f1d): MCISamplingFunctionInterface(1, 1) {_f1d=f1d;}

  void samplingFunction(const double * in, double * protovalue){
    protovalue[0] = 1.; // set something reasonable here
  }

  double getAcceptance(){
    return this->getProtoNew(0) / this->getProtoOld(0);
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

      _ffnn->setInput(1, in);
      _ffnn->FFPropagate();
      double fx = _ftarget->f(in[0]);
      out[0] = pow(_ffnn->getOutput(1) - fx, 2);

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

    _ffnn->setInput(1, in);
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

    //initial walker position and step size
    _x0 = 0.5*(irange[0] + irange[1]);
    _step0 = 0.05*(irange[1] - irange[0]);

    //setup MCIs for distance function and gradient
    _mcif = new MCI(1);
    _mcig = new MCI(1);

    MCIObservableFunctionInterface * fobs = new NNFitDistance1D(_ffnn, _ftarget);
    MCIObservableFunctionInterface * gobs = new NNFitGradient1D(_ffnn, _ftarget);
    MCISamplingFunctionInterface * fts = new SamplingFunction1D(_ftarget);

    _mcif->addSamplingFunction(fts);
    _mcig->addSamplingFunction(fts);

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

// creates and holds the necessary data copies for multiple fit threads in parallel
class NNFitter1D{
private:
  int _nhlayer;
  int * _nhunits;
  long _nmc;
  double * _irange;

  Function1D * _ftarget;
  FeedForwardNeuralNetwork * _ffnn;
  FitNN1D * _fitnn;
  ConjGrad * _conj;

public:
  NNFitter1D(Function1D * ftarget, int &nhlayer, int * nhunits, long &nmc, double * irange){
    _ftarget = ftarget;
    _nhlayer = nhlayer;
    _nhunits = nhunits;
    _nmc = nmc;
    _irange = irange;
  }

  void init() { // capsule this away for threading
    _ffnn = new FeedForwardNeuralNetwork(2, _nhunits[0], 2);
    for (int i = 1; i<_nhlayer; ++i) _ffnn->pushHiddenLayer(_nhunits[i]);
    _ffnn->connectFFNN();
    _ffnn->addVariationalFirstDerivativeSubstrate();

    _fitnn = new FitNN1D(_ffnn, _ftarget, _nmc, _irange);

    _conj = new ConjGrad(_fitnn);
    for(int i = 0; i<_ffnn->getNBeta(); ++i) _conj->setX(i, _ffnn->getBeta(i));
  }

  void findFit() {
    _conj->findMin();
  }

  double getFitDistance(){
    return 0.;
  }

  FeedForwardNeuralNetwork * getFFNN() {return _ffnn;}
  FitNN1D * getFitNN() {return _fitnn;}
  ConjGrad * getConj() {return _conj;}
};

// code which runs in parallel
void *findFit_thread(void * voidPtr) {
  NNFitter1D * fitter = static_cast<NNFitter1D*>(voidPtr);
  printf("This is findFit_thread()\n");
  fitter->init();
  fitter->findFit();
  printf("Still here\n");
  pthread_exit(NULL);
}


int main() {
  using namespace std;

  int nl, nhl, nhu[2];
  long nmc = 40000l;
  //double epsx = 0.01;
  double irange[2];
  irange[0] = -10.;
  irange[1] = 10.;

  cout << "Let's start by creating a Feed Forward Artificial Neural Netowrk (FFANN)" << endl;
  cout << "========================================================================" << endl;
  cin.ignore();

  cout << "How many units should the first hidden layer(s) have? ";
  cin >> nhu[0];
  cout << "How many units should the second hidden layer(s) have? (<=1 for no second hidden layer)";
  cin >> nhu[1];

  nl = (nhu[1]>1)? 4:3;
  nhl = nl-2;
  cout << "We generate a FFANN with " << nl << " layers and 2, " << nhu[0];
  if (nhu[1]>1) { cout << ", " << nhu[1];}
  cout << ", 2 units respectively" << endl;
  cout << "========================================================================" << endl << endl;
  cin.ignore();

  /*
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
  */

  // NON I/O CODE
  int nthread = 4;
  pthread_t * my_threads = new pthread_t[nthread];
  int ret[nthread];

  Gaussian * gauss = new Gaussian(0.5,0.);
  NNFitter1D * fit_list[nthread];
  //

  cout << "Now let CG find the optimal betas..." << endl << endl;
  cin.ignore();

  // NON I/O CODE
  // Preparing list of FitStructs for threads
  for (int i=0; i<nthread; ++i){
    fit_list[i] = new NNFitter1D(gauss, nhl, nhu, nmc, irange);
  }
  printf("In main: creating threads\n");
  for (int i=0; i<nthread; ++i){
    ret[i] =  pthread_create(&my_threads[i], NULL, &findFit_thread, fit_list[i]);
    if(ret[i] != 0) {
      printf("Error: pthread_create() failed\n");
      exit(EXIT_FAILURE);
    }
  }
  for (int i=0; i<nthread; ++i) pthread_join(my_threads[i], (void **) &ret[i]);
  printf("Threads joined\n");
  //

  cout << "Done." << endl;
  cout << "Now copy over the optimal betas from CG to the NN and evaluate distance." << endl << endl;
  cin.ignore();

  // NON I/O CODE
  for (int i=0; i<nthread; ++i) {
    for(int j = 0; j<fit_list[i]->getFFNN()->getNBeta(); ++j) {
      fit_list[i]->getFFNN()->setBeta(j, fit_list[i]->getConj()->getX(j));
    }
  }
  //

  cout << "Finally the NN is fitted!" << endl << endl;
  cin.ignore();

  cout << "Let us compare the NN to target function values:" << endl;
  cin.ignore();

  double x=-10., dx=1.;
  for(int i=0; i<21; ++i) {
    fit_list[0]->getFFNN()->setInput(1, &x);
    fit_list[0]->getFFNN()->FFPropagate();
    cout << "x: " << x << " f(x): " << gauss->f(x) << " nn(x): " << fit_list[0]->getFFNN()->getOutput(1) << endl;
    x+=dx;
  }
  cout << endl;
  /*
  for (int i=0; i<nthread; ++i) {
    delete conj_list[i];
    delete fitnn_list[i];
    delete ffnn_list[i];
    }*/
  //delete conj_list;
  //delete fitnn;
  //delete gauss;
  //delete ffnn;

  // end
  return 0;
}
