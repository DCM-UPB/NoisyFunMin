#include <iostream>
#include <cmath>
#include <math.h>
#include <pthread.h>

#include "NoisyFunction.hpp"
#include "ConjGrad.hpp"
#include "MCIntegrator.hpp"
#include "FeedForwardNeuralNetwork.hpp"
#include "PrintUtilities.hpp"

/*
Example: Fit FFNN via NFM's Conjugate Gradient

In this example we want to demonstrate how a FFNN can be trained via Conjuage Gradient (CG) to fit a target function (in this case a gaussian).

This is done by minimizing the quadratic distance of NN to target function via class CG from NFM. The distance and gradient functions are integrated by MCI.
Hence, we have to create the NN and target function, MCIObservableFunctionInterfaces for distance and gradient, and put everything together into a NoisyFunctionWithGradient for CG.

The whole process is handled by the NNFitter1D class. As a little speciality, the fitting is done multiple times in parallel threads and the best fit overall is presented as result.
*/


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

      _ffnn->setInput(1, in);
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

    _ffnn->setInput(1, in);
    _ffnn->FFPropagate();

    double nnout = _ffnn->getOutput(1);
    double fx = _ftarget->f(in[0]);

    for(int i=0; i<_nobs; ++i){
      out[i] = 2.*(nnout - fx)*_ffnn->getVariationalFirstDerivative(1, i);
    }
  }
};

// NFwG for minimization of distance between NN and target function
class FitNN1D: public NoisyFunctionWithGradient{
private:
  FeedForwardNeuralNetwork * _ffnn;
  Function1D * _ftarget;
  long _nmc;

  double ** _intrange;
  MCI * _mcif, * _mcig;
  MCIObservableFunctionInterface * _fobs, * _gobs;


public:
  FitNN1D(FeedForwardNeuralNetwork * ffnn, Function1D * ftarget, const long &nmc, double * irange): NoisyFunctionWithGradient(ffnn->getNBeta()){
    _ffnn = ffnn;
    _ftarget = ftarget;
    _nmc = nmc;

    //setup MCIs for distance function and gradient
    _mcif = new MCI(1);
    _mcig = new MCI(1);

    _fobs = new NNFitDistance1D(_ffnn, _ftarget);
    _gobs = new NNFitGradient1D(_ffnn, _ftarget);

    _mcif->addObservable(_fobs);
    _mcig->addObservable(_gobs);

    _intrange = new double*[1];
    _intrange[0] = irange;
    _mcif->setIRange(_intrange);
    _mcig->setIRange(_intrange);

  }

  ~FitNN1D() {
    delete _mcif;
    delete _mcig;
    delete _fobs;
    delete _gobs;
    delete _intrange;
  }


  void f(const double * in, double &f, double &df){

    //set new NN betas
    for (int i=0; i<_ffnn->getNBeta(); ++i){
      _ffnn->setBeta(i, in[i]);
    }

    // integrate
    _mcif->integrate(_nmc, &f, &df);
  }

  void grad(const double * in, double * g, double * dg){
    int i, nbeta = _ffnn->getNBeta();

    //set new NN betas
    for (i=0; i<nbeta; ++i){
      _ffnn->setBeta(i, in[i]);
    }

    // integrate
    _mcig->integrate(_nmc, g, dg);
  }

};

// creates instances and holds the necessary data for multiple fit threads in parallel
class NNFitter1D {
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
  NNFitter1D(Function1D * ftarget, const int &nhlayer, int * nhunits, const long &nmc, double * irange) {
    _ftarget = ftarget;
    _nhlayer = nhlayer;
    _nhunits = nhunits;
    _nmc = nmc;
    _irange = irange;
  }

  ~NNFitter1D(){
    delete _ffnn;
    delete _fitnn;
    delete _conj;
  }

  void init() { // capsule this init part away for threading
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

  // compute fit distance for CG's best betas
  double getFitDistance() {
    double betas[_ffnn->getNBeta()];
    double f,df;

    for (int i=0; i<_ffnn->getNBeta(); ++i) {
      betas[i] = _conj->getX(i);
    }
    _fitnn->f(betas, f, df);
    return f;
  }

  // compare NN to target at nx points starting from x=x0 in increments dx
  void compareFit(const double x0, const double dx, const int nx) {
    using namespace std;
    double x=x0;
    for(int i=0; i<nx; ++i) {
      _ffnn->setInput(1, &x);
      _ffnn->FFPropagate();
      cout << "x: " << x << " f(x): " << _ftarget->f(x) << " nn(x): " << _ffnn->getOutput(1) << endl;
      x+=dx;
    }
    cout << endl;
  }

  // print output of fitted NN to file
  void printFitOutput(int input_i, int output_i, const double &min, const double &max, const int &npoints) {
    double * base_input = new double[_ffnn->getNInput()];
    writePlotFile(_ffnn, base_input, input_i, output_i, min, max, npoints, "getOutput", "v.txt");
    delete [] base_input;
  }

  // store fitted NN in file
  void printFitNN() {_ffnn->storeOnFile("nn.txt");}

  FeedForwardNeuralNetwork * getFFNN() {return _ffnn;}
  FitNN1D * getFitNN() {return _fitnn;}
  ConjGrad * getConj() {return _conj;}
};

// code which runs in parallel
void findFit(void * voidPtr) {
  using namespace std;
  NNFitter1D * fitter = static_cast<NNFitter1D*>(voidPtr);
  fitter->init();
  try {
    fitter->findFit();
  }
  catch (runtime_error e) {
    cout << "Warning: Fit thread aborted because of too long bracketing." << endl;
  }
}

void *findFit_thread(void * voidPtr) {
  findFit(voidPtr);
  pthread_exit(NULL);
}


int main() {
  using namespace std;

  int nl, nhl, nhu[2], nthread;
  const long nmc = 40000l;
  double irange[2];
  irange[0] = -10.;
  irange[1] = 10.;

  cout << "Let's start by creating a Feed Forward Artificial Neural Network (FFANN)" << endl;
  cout << "========================================================================" << endl;
  cout << endl;
  cout << "How many units should the first hidden layer(s) have? ";
  cin >> nhu[0];
  cout << "How many units should the second hidden layer(s) have? (<=1 for none) ";
  cin >> nhu[1];
  cout << endl;

  // NON I/O CODE
  nl = (nhu[1]>1)? 4:3;
  nhl = nl-2;
  //

  cout << "We generate a FFANN with " << nl << " layers and 2, " << nhu[0];
  if (nhu[1]>1) { cout << ", " << nhu[1];}
  cout << ", 2 units respectively" << endl;
  cout << "========================================================================" << endl;
  cout << endl;
  cout << "In the following we use CG to minimize the mean-squared-distance of NN vs. target function, i.e. find optimal betas." << endl;
  cout << endl;
  cout << "How many fitting threads do you want to spawn? ";
  cin >> nthread;
  cout << endl << endl;
  cout << "Now we spawn fitting threads and find the best fit ... " << endl;;

  // NON I/O CODE

  // Preparing list of NNFitters for the threads
  pthread_t * my_threads = new pthread_t[nthread];
  int ret[nthread];

  Gaussian * gauss = new Gaussian(0.5,0.);

  NNFitter1D * fit_list[nthread];
  for (int i=0; i<nthread; ++i){
    fit_list[i] = new NNFitter1D(gauss, nhl, nhu, nmc, irange);
  }

  // Spawn minimzation threads
  for (int i=0; i<nthread; ++i) {
    ret[i] =  pthread_create(&my_threads[i], NULL, &findFit_thread, fit_list[i]);
    if(ret[i] != 0) {
      printf("Error: pthread_create() failed\n");
      exit(EXIT_FAILURE);
    }
  }

  // Join the threads
  for (int i=0; i<nthread; ++i) pthread_join(my_threads[i], (void **) &ret[i]);

  // Find best fit index
  int bfi = 0;
  double fdist, bdist = fit_list[0]->getFitDistance();

  for(int i=1; i<nthread; ++i) {
    fdist = fit_list[i]->getFitDistance();
    if(fdist < bdist) {
      bfi = i;
      bdist = fdist;
    }
  }

  //

  cout << "Done." << endl;
  cout << "========================================================================" << endl;
  cout << endl;
  cout << "Finally we compare the best fit NN to the target function:" << endl << endl;

  // NON I/O CODE
  fit_list[bfi]->compareFit(-10., 1., 21);
  //

  cout << endl;
  cout << "And print the output/NN to a file. The end." << endl;

  // NON I/O CODE
  fit_list[bfi]->printFitOutput(0, 1, -10, 10, 200);
  fit_list[bfi]->printFitNN();
  //

  // cleanup

  delete [] my_threads;
  delete gauss;
  for (int i=0; i<nthread; ++i) delete fit_list[i];

  // end
  return 0;
}

/* no-thread main for debug
int main() {
  using namespace std;

  int nl, nhl, nhu[2];
  nhu[0] = 15;
  nhu[1] = 0;
  nl = (nhu[1]>1)? 4:3;
  nhl = nl-2;


  const long nmc = 1000l;
  double irange[2];
  irange[0] = -10.;
  irange[1] = 10.;

  Gaussian * gauss = new Gaussian(0.5,0.);

  NNFitter1D * fitter = new NNFitter1D(gauss, nhl, nhu, nmc, irange);
  fitter->init();
  fitter->findFit();
  fitter->printFitOutput(0, 1, -10, 10, 200);



  delete gauss;
  delete fitter;
};
*/
