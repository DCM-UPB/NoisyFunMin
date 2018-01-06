#include "ConjGrad.hpp"

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

void ConjGrad::writeCurrentXInLog(){
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
   
   delete log_manager;
}


void ConjGrad::writeDirectionInLog(const double * direction, const double * directionerror){
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
   
   delete log_manager;
}


void ConjGrad::reportMeaninglessGradientInLog(){
   using namespace std;
   
   NFMLogManager * log_manager = new NFMLogManager();
      
   stringstream s;
   s << endl << "gradient seems to be meaningless, i.e. its error is too large" << endl;
   s << flush;
   log_manager->writeOnLog(s.str());   
   
   delete log_manager;
}



// --- Minimization

void ConjGrad::findMin()
{
   using namespace std;
   
   NFMLogManager * log_manager = new NFMLogManager();
   log_manager->writeOnLog("\nBegin ConjGrad::findMin() procedure\n");

   //check if the gradient has been provided before starting the minimization
   if ( !_flaggradtargetfun )
   {
      cout << "ERROR ConjGrad.findMin() : The gradient is required for this method, but it was not provided" << endl << endl;
      exit(-1);
   }

   ////check if the starting X is inside the domain
   //if ( _flagindomain )
   //{
   //   if ( !_indomain(_ndim,_x) )
   //   {
   //      cout << "ERROR ConjGrad.findMin() : Starting X out of the domain" << std::endl << std::endl;
   //      exit(-2);
   //   }
   //}

   //initialize the gradients
   double * gradold = new double[_ndim];
   double * gradnew = new double[_ndim];
   double * graderr = new double[_ndim];
   double * conjvold = new double[_ndim];
   double * conjvnew = new double[_ndim];
   
   this->_gradtargetfun->grad(_x->getX(), gradold, graderr);
   
   this->writeCurrentXInLog();
   
   if (this->_meaningfulGradient(gradold, graderr))
   {
      int i;
      //interested in following -gradient
      for (i=0; i<_ndim; ++i){ gradold[i]=-gradold[i]; }
      //inizialize the conjugate vectors
      for (i=0; i<_ndim; ++i){ conjvold[i]=gradold[i]; }
      this->writeDirectionInLog(conjvold, graderr);
      //find new position
      double deltatargetfun, deltax;
      this->findNextX(conjvold,deltatargetfun,deltax);
      
      this->writeCurrentXInLog();

      //begin the minimization loop
      double scalprodold, scalprodnew, ratio;
      //cout << "deltatargetfunction = " << deltatargetfun << "   " << _epstargetfun << endl;
      //cout << "deltax = " << deltax << "   " << _epsx << endl << endl;
      while ( ( deltatargetfun>=_epstargetfun ) && (deltax>=_epsx) )
      {
         //cout << "x is in " << getX(0) << "   " << getX(1) << "   " << getX(2) << endl << endl;
         //evaluate the new gradient
         this->_gradtargetfun->grad(_x->getX(),gradnew,graderr);
         for (i=0; i<_ndim; ++i){ gradnew[i]=-gradnew[i]; }
         if (!this->_meaningfulGradient(gradnew, graderr))
         {
            this->reportMeaninglessGradientInLog();
            break;
         }
         // compute the direction to follow for finding the next x
         // if _use_conjgrad == true -> Conjugate Gradient
         // else -> Steepest Descent
         if (_use_conjgrad){
            //determine the new conjugate vector
            scalprodnew=0.;
            for (i=0; i<_ndim; ++i){ scalprodnew+=gradnew[i]*gradnew[i]; }
            scalprodold=0.;
            for (i=0; i<_ndim; ++i){ scalprodold+=gradold[i]*gradold[i]; }
            ratio=scalprodnew/scalprodold;
            for (i=0; i<_ndim; ++i){ conjvnew[i]=gradnew[i]+conjvold[i]*ratio; }
         } else {
            // simply use as conjugate gradient the gradient (i.e. make a steepest descent!)
            for (i=0; i<_ndim; ++i){ conjvnew[i]=gradnew[i]; }
         }
         this->writeDirectionInLog(conjvnew, graderr);
         //find new position
         this->findNextX(conjvnew,deltatargetfun,deltax);
         //cout << "deltatargetfunction = " << deltatargetfun << "   " << _epstargetfun << endl;
         //cout << "deltax = " << deltax << "   " << _epsx << endl << endl;
         
         this->writeCurrentXInLog();
      }
   }
   
   log_manager->writeOnLog("\nEnd ConjGrad::findMin() procedure\n"); 

   //free memory
   delete [] conjvnew;
   delete [] conjvold;
   delete [] gradnew;
   delete [] gradold;
   delete [] graderr;
   
   delete log_manager;
}


// --- Internal methods

void ConjGrad::findNextX(const double * dir, double &deltatargetfun, double &deltax)
{
   using namespace std;
   
   //project the original multidimensional wave function into a one-dimensional function
   FunProjection1D * proj1d = new FunProjection1D(this->_ndim, this->_x->getX(), dir, this->_targetfun);
   //determine the initial bracket
   NoisyFunctionValue a(1), b(1), c(1);
   a.setX(0.);
   double newf, dnewf;
   proj1d->f(a.getX(), newf, dnewf);
   a.setF(newf, dnewf);
   nfm::findBracket(proj1d, a, b, c);
   //find the minimum in the bracket
   nfm::parabgoldMinimization(proj1d, this->_epstargetfun, a, b, c);
   //get the x corresponding to the found b
   proj1d->getVecFromX(b.getX(0),_x->getX());
   _x->setF(b.getF(), b.getDf());
   //compute the two deltas
   deltatargetfun=std::abs(b.getF()-newf)-dnewf-b.getDf();
   deltax=0.;
   for (int i=0; i<this->_ndim; ++i)
   {
      deltax+=dir[i]*dir[i];
   }
   deltax=std::abs(b.getX(0)*sqrt(deltax)); //*NORMA VETTORE dir
   //cout << "x is in " << this->getX(0) << "   " << this->getX(1) << "   " << this->getX(2) << endl;
   delete proj1d;
}


