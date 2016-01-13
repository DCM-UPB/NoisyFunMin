#ifndef FUN_PROJECTION_1D
#define FUN_PROJECTION_1D

#include "NoisyFunction.hpp"


class FunProjection1D: public NoisyFunction
{
protected:
   double * _p0;   //starting point
   double * _direction;   //direction
   double * _vec;  //vector used internally
   NoisyFunction * _mdf;  //multidimensional function that must be projected
   int _originalndim;  //dimension of the multidimensional function
   
public:
   FunProjection1D(const int & originalndim, const double * p0, const double * direction, NoisyFunction * mdf): NoisyFunction(1)
   {
      _originalndim=originalndim;
      _p0=new double[_originalndim];
      for (int i=0; i<_originalndim; ++i)
      {
         _p0[i]=p0[i];
      }
      _direction=new double[_originalndim];
      for (int i=0; i<_originalndim; ++i)
      {
         _direction[i]=direction[i];
      }
      _vec=new double[_originalndim];
      _mdf=mdf;
   }
   ~FunProjection1D();
   
   void f(const double *x, double &f, double &df);  //projected one-dimensional function

   void getVecFromX(const double &x, double *vec);
   
};

#endif
