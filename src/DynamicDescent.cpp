#include "DynamicDescent.hpp"

#include "LogNFM.hpp"
#include "FunProjection1D.hpp"
#include "NoisyFunctionValue.hpp"
#include "1DTools.hpp"

#include "stdlib.h"
#include <iostream>
#include <sstream>
#include <cmath>
#include <stdexcept>




// --- Log

void DynamicDescent::writeCurrentXInLog(){
   using namespace std;
   
   NFMLogManager * log_manager = new NFMLogManager();
      
   stringstream s;
   s << endl << "x:\n";
   for (int i=0; i<_x->getNDim(); ++i){
      s << _x->getX(i) << "    ";
   }
   s << endl << "    ->    value = " << _x->getF() << " +- " << _x->getDf() << endl;
   s << flush;
      
   log_manager->writeOnLog(s.str());
}


void DynamicDescent::writeDirectionInLog(const double * direction, const double * directionerror){
   using namespace std;
   
   NFMLogManager * log_manager = new NFMLogManager();
      
   stringstream s;
   s << endl << "direction (and error):\n";
   for (int i=0; i<_x->getNDim(); ++i){
      s << direction[i] << " (" << directionerror[i] << ")    ";
   }
   s << endl;
   s << flush;
   log_manager->writeOnLog(s.str());
}


void DynamicDescent::writeInertiaInLog(){
   using namespace std;
   
   NFMLogManager * log_manager = new NFMLogManager();
      
   stringstream s;
   s << endl << "inertia: " << _inertia << endl;
   s << flush;
   log_manager->writeOnLog(s.str());
}


void DynamicDescent::reportMeaninglessGradientInLog(){
   using namespace std;
   
   NFMLogManager * log_manager = new NFMLogManager();
      
   stringstream s;
   s << endl << "gradient seems to be meaningless, i.e. its error is too large" << endl;
   s << flush;
   log_manager->writeOnLog(s.str());   
}



// --- Minimization

void DynamicDescent::findMin()
{
   using namespace std;
   
   NFMLogManager * log_manager = new NFMLogManager();
   log_manager->writeOnLog("\nBegin DynamicDescent::findMin() procedure\n");

   //check if the gradient has been provided before starting the minimization
   if ( !_flaggradtargetfun )
   {
      cout << "ERROR DynamicDescent.findMin() : The gradient is required for this method, but it was not provided" << endl << endl;
      exit(-1);
   }

   //initialize the gradients
   double * gradold = new double[_ndim];
   double * gradnew = new double[_ndim];
   double * graderr = new double[_ndim];
   
   this->_gradtargetfun->grad(_x->getX(), gradold, graderr);
   
   _inertia = 0.;
   for (int i=0; i<this->_ndim; ++i)
   {
      _inertia+=gradold[i]*gradold[i];
   }
   _inertia = 1. / sqrt(_inertia);
   
   this->writeCurrentXInLog();
   
   if (this->_meaningfulGradient(gradold, graderr))
   {
      int i;
      //interested in following -gradient
      for (i=0; i<_ndim; ++i){ gradold[i]=-gradold[i]; }
      this->writeDirectionInLog(gradold, graderr);
      //find new position
      double deltatargetfun, deltax;
      this->findNextX(gradold, deltatargetfun, deltax);
      this->writeCurrentXInLog();
      //begin the minimization loop
      while ( ( deltatargetfun>=_epstargetfun ) && (deltax>=_epsx) )
      {
         //evaluate the new gradient
         this->_gradtargetfun->grad(_x->getX(),gradnew,graderr);
         for (i=0; i<_ndim; ++i){ gradnew[i]=-gradnew[i]; }
         if (!this->_meaningfulGradient(gradnew, graderr))
         {
            this->reportMeaninglessGradientInLog();
            break;
         }
         this->writeDirectionInLog(gradnew, graderr);
         this->findNextX(gradnew,deltatargetfun,deltax);
         this->writeCurrentXInLog();
      }
   }
   
   log_manager->writeOnLog("\nEnd DynamicDescent::findMin() procedure\n"); 

   //free memory
   delete [] gradnew;
   delete [] gradold;
   delete [] graderr;
   
   delete log_manager;
}


// --- Internal methods

void DynamicDescent::findNextX(const double * dir, double &deltatargetfun, double &deltax)
{
   using namespace std;
   
   // compute the normalized direction vector
   double * norm_dir = new double[_ndim];
   double sum = 0.;
   for (int i=0; i<_ndim; ++i){
      sum += dir[i]*dir[i];
   }
   sum = sqrt(sum);
   for (int i=0; i<_ndim; ++i){
      norm_dir[i] = dir[i]/sum;
   }
   // compute the dot product between the normalized direction and the old normalized direction
   double old_new_direction_dot_product = 0;
   for (int i=0; i<_ndim; ++i){
      old_new_direction_dot_product += _old_norm_direction[i] * norm_dir[i];
   }
   // update the inertia and write it in the log
   _inertia = _inertia + 0.5 * _inertia * old_new_direction_dot_product;
   this->writeInertiaInLog();
   
   // update _x
   for (int i=0; i<_ndim; ++i){
      _x->setX(i, _x->getX(i) + _inertia * dir[i]);
   }
   
   // compute deltax
   deltax=0.;
   for (int i=0; i<this->_ndim; ++i)
   {
      deltax+=dir[i]*dir[i];
   }
   deltax=_inertia * sqrt(deltax);
   // compute the new value of the target function
   double newf, newdf;
   _gradtargetfun->f(_x->getX(), newf, newdf);
   // compute the deltatargetfunction
   deltatargetfun=std::abs(_x->getF()-newf)-newdf-_x->getDf();
   // store the new function value in _x
   _x->setF(newf, newdf);
   
   // store the direction for the next iteration
   for (int i=0; i<this->_ndim; ++i){
      _old_norm_direction[i] = norm_dir[i];
   }
   
   delete[] norm_dir;
}


