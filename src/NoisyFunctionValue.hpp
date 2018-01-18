#ifndef NOISY_FUNCTION_VALUE
#define NOISY_FUNCTION_VALUE

// function: R^(ndim) -> R   with noise (i.e. the value of the function has an associated error)
// this class implements one of its value (the input and its output, i.e. x and f(x) with associated error)
class NoisyFunctionValue
{
   protected:
      int _ndim;  //dimensionality
      double * _x;   //point where the function has been evaluated
      double _f;   //value of the function
      double _df;  //error associated to f

   public:
      NoisyFunctionValue(int ndim);
      ~NoisyFunctionValue();

      //Setters
      void setX(const double &x){_x[0]=x;}
      void setX(const double * x);
      void setX(const int &i, const double &x){_x[i]=x;}
      void setF(const double &f, const double &df);

      //Getters
      int getNDim(){return _ndim;}
      double getX(const int &idx){return _x[idx];}
      double * getX(){return _x;}
      double getF(){return _f;}
      double getDf(){return _df;}

      //Operation
      bool operator< (const NoisyFunctionValue &val);
      bool operator>= (const NoisyFunctionValue &val){return !((*this)<val);}
      bool operator> (const NoisyFunctionValue &val);
      bool operator<= (const NoisyFunctionValue &val){return !((*this)>val);}
      bool operator== (const NoisyFunctionValue &val);
      NoisyFunctionValue& operator= (const NoisyFunctionValue &val);
};


#endif
