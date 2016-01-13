#include "ConjGrad.hpp"

#include "FunProjection1D.hpp"
#include "NoisyFunctionValue.hpp"
#include "1DTools.hpp"

#include "stdlib.h"
#include <iostream>
#include <cmath>

// --- Minimization

void ConjGrad::findMin()
{
   using namespace std;

   //check if the gradient has been provided before starting the minimization
   if ( !_flaggradtargetfun )
   {
      cout << "ERROR ConjGrad.findMin() : The gradient is required for this method, but it was not provided" << std::endl << std::endl;
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
   this->_gradtargetfun->grad(_x->getX(),gradold, graderr);
   if (this->_meaningfulGradient(gradold, graderr))
   {
      int i;
      //interested in following -gradient
      for (i=0; i<_ndim; ++i){ gradold[i]=-gradold[i]; }
      //inizialize the conjugate vectors
      for (i=0; i<_ndim; ++i){ conjvold[i]=gradold[i]; }
      //find new position
      double deltatargetfun, deltax;
      this->findNextX(conjvold,deltatargetfun,deltax);

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
            //cout << "Gradient meaningless!" << endl;
            //cout << gradnew[0] << "   " << gradnew[1] << "   " << gradnew[2] << endl;
            //cout << graderr[0] << "   " << graderr[1] << "   " << graderr[2] << endl;
            break;
         }
         //determine the new conjugate vector
         scalprodnew=0.;
         for (i=0; i<_ndim; ++i){ scalprodnew+=gradnew[i]*gradnew[i]; }
         scalprodold=0.;
         for (i=0; i<_ndim; ++i){ scalprodold+=gradold[i]*gradold[i]; }
         ratio=scalprodnew/scalprodold;
         for (i=0; i<_ndim; ++i){ conjvnew[i]=gradnew[i]+conjvold[i]*ratio; }
         //find new position
         this->findNextX(conjvnew,deltatargetfun,deltax);
         //cout << "deltatargetfunction = " << deltatargetfun << "   " << _epstargetfun << endl;
         //cout << "deltax = " << deltax << "   " << _epsx << endl << endl;
      }
   }

   //free memory
   delete [] conjvnew;
   delete [] conjvold;
   delete [] gradnew;
   delete [] gradold;
   delete [] graderr;
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
   proj1d->f(a.getX(),newf,dnewf);
   a.setF(newf,dnewf);
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


