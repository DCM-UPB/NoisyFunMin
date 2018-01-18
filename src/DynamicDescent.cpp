#include "DynamicDescent.hpp"

#include "LogNFM.hpp"
#include "FunProjection1D.hpp"
#include "1DTools.hpp"

#include "stdlib.h"
#include <iostream>
#include <sstream>
#include <cmath>
#include <stdexcept>




// --- Log

void DynamicDescent::writeCurrentXInLog(){
   using namespace std;
   
   NFMLogManager log_manager = NFMLogManager();
      
   stringstream s;
   s << endl << "x:\n";
   for (int i=0; i<_x->getNDim(); ++i){
      s << _x->getX(i) << "    ";
   }
   s << endl << "    ->    value = " << _x->getF() << " +- " << _x->getDf() << endl;
   s << flush;
   log_manager.writeOnLog(s.str());   
}


void DynamicDescent::writeDirectionInLog(const double * grad, const double * graderror){
   using namespace std;
   
   NFMLogManager log_manager = NFMLogManager();
      
   stringstream s;
   s << endl << "grad (and error):\n";
   for (int i=0; i<_x->getNDim(); ++i){
      s << -grad[i] << " (" << graderror[i] << ")    ";
   }
   s << endl;
   s << flush;
   log_manager.writeOnLog(s.str());
}


void DynamicDescent::writeInertiaInLog(){
   using namespace std;
   
   NFMLogManager log_manager = NFMLogManager();
      
   stringstream s;
   s << endl << "inertia: " << _inertia << endl;
   s << flush;
   log_manager.writeOnLog(s.str());
}


void DynamicDescent::reportMeaninglessGradientInLog(){
   using namespace std;
   
   NFMLogManager log_manager = NFMLogManager();
      
   stringstream s;
   s << endl << "gradient seems to be meaningless, i.e. its error is too large" << endl;
   s << flush;
   log_manager.writeOnLog(s.str());
}


void DynamicDescent::writeOldValuesInLog(){
   using namespace std;
   
   NFMLogManager log_manager = NFMLogManager();
      
   stringstream s;
   s << endl << "last values:    ";
   for (list<NoisyFunctionValue *>::iterator it=_old_values.begin(); it!=_old_values.end(); ++it){
      s << (*it)->getF() << " +- " << (*it)->getDf() << "    ";
   }
   s << endl;
   s << "equal to first element? ";
   for (list<NoisyFunctionValue *>::iterator it=_old_values.begin(); it!=_old_values.end(); ++it){
      if (it != _old_values.begin()){
         s << ((**it) == (**_old_values.begin())) << "    ";
      }
      
   }
   s << endl;
   log_manager.writeOnLog(s.str());
}



// --- Minimization

void DynamicDescent::findMin(){
   using namespace std;
      
   NFMLogManager * log_manager = new NFMLogManager();
   log_manager->writeOnLog("\nBegin DynamicDescent::findMin() procedure\n");
   
   //initialize the gradients
   double * grad = new double[_ndim];
   double * graderr = new double[_ndim];
   
   // compute the current value
   double newf, newdf;
   this->_gradtargetfun->f(_x->getX(), newf, newdf);
   _x->setF(newf, newdf);
   this->writeCurrentXInLog();
   
   
   //begin the minimization loop
   int cont = 0;
   while ( shouldContinueDescent() )
   {
      // compute the gradient
      this->_gradtargetfun->grad(_x->getX(), grad, graderr);
      this->writeDirectionInLog(grad, graderr);
      
      // if it is the first iteration, initialise the inertia
      if (cont == 0){
         _inertia = 0.;
         for (int i=0; i<this->_ndim; ++i)
         {
            _inertia+=grad[i]*grad[i];
         }
         _inertia = _ndim / sqrt(_inertia);
      }
      
      // find the next position
      this->findNextX(grad);
      this->writeCurrentXInLog();
      
      cont ++;
   }
   
   log_manager->writeOnLog("\nEnd DynamicDescent::findMin() procedure\n"); 
   
   //free memory
   delete [] grad;
   delete [] graderr;
   
   //delete log_manager;
}


// --- Internal methods

void DynamicDescent::findNextX(const double * grad)
{
   using namespace std;
   
   // compute the normalized gradection vector
   double * norm_grad = new double[_ndim];
   double sum = 0.;
   for (int i=0; i<_ndim; ++i){
      sum += grad[i]*grad[i];
   }
   sum = sqrt(sum);
   for (int i=0; i<_ndim; ++i){
      norm_grad[i] = grad[i]/sum;
   }
   // compute the dot product between the normalized gradection and the old normalized gradection
   double old_new_direction_dot_product = 0;
   for (int i=0; i<_ndim; ++i){
      old_new_direction_dot_product += _old_norm_direction[i] * norm_grad[i];
   }
   // update the inertia and report it in the log
   _inertia = _inertia + 0.5 * _inertia * old_new_direction_dot_product;
   this->writeInertiaInLog();
   
   // update _x
   for (int i=0; i<_ndim; ++i){
      _x->setX(i, _x->getX(i) - _inertia * grad[i]);
   }
   
   // compute the new value of the target function
   double newf, newdf;
   cout << "A" << endl;
   _gradtargetfun->f(_x->getX(), newf, newdf);
   cout << "B" << endl;

   // store the new function value in _x
   _x->setF(newf, newdf);
   
   // store the grad for the next iteration
   for (int i=0; i<this->_ndim; ++i){
      _old_norm_direction[i] = norm_grad[i];
   }
   
   delete[] norm_grad;
}


bool DynamicDescent::shouldContinueDescent(){
   using namespace std;
      
   NoisyFunctionValue * v = new NoisyFunctionValue(_x->getNDim());
   *v = *_x;
   _old_values.push_front(v);
      
   if (_old_values.size() > N_CONSTANT_VALUES_CONDITION_FOR_STOP){
      _old_values.pop_back();
   }
   
   writeOldValuesInLog();
   
   if (_old_values.size() == N_CONSTANT_VALUES_CONDITION_FOR_STOP){
      
      for (list<NoisyFunctionValue *>::iterator it = _old_values.begin(); it != _old_values.end(); ++it){
         if (it != _old_values.begin()){
            if (! (**it == **_old_values.begin()) ){
               return true;
            }
         }
      }
      
      NFMLogManager log_manager = NFMLogManager();
      log_manager.writeOnLog("\nValues have stabilised, interrupting minimisation procedure.\n");
            
      return false;
   }
      
   return true;
}


