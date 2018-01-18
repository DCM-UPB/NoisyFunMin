#include "NoisyFunMin.hpp"

#include <cmath>
#include <iostream>



// --- Protected methods

bool NFM::_meaningfulGradient(const double * grad, const double * graderr)
{
   for (int i=0; i<_ndim; ++i)
   {
      if (std::abs(grad[i])>graderr[i]) return true;
   }
   return false;
}


// --- Getters




// --- Setters

//void NFM::setDomain(nfm::DomainFun domain)
//{
//   _indomain = domain;
//   _flagindomain = true;
//}


void NFM::setGradientTargetFun(NoisyFunctionWithGradient * grad)
{
   _gradtargetfun = grad;
   _flaggradtargetfun = true;
}


void NFM::setX(const double * x)
{
   _x->setX(x);
}


void NFM::setX(const int &i, const double &x)
{
   _x->setX(i, x);
}


// --- Constructor and destructor

NFM::NFM(NoisyFunction * targetfun)
{
   //set ndim and the target function
   _targetfun = targetfun;
   _ndim = _targetfun->getNDim();
   //allocate and initialize x
   _x = new NoisyFunctionValue(_ndim);
   for (int i=0; i<_ndim; ++i){ _x->setX(i,0.); }
   //gradient of the target function
   _gradtargetfun = 0;
   _flaggradtargetfun = false;
   //optimization's domain
   //_indomain = 0;
   //_flagindomain = false;
   //eps
   _epstargetfun=0.;
   _epsx=0.;
}


NFM::~NFM()
{
   //deallocate x
   delete _x;
}
