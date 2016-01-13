#include "1DTools.hpp"

#include <iostream>

namespace nfm
{
   
   void findBracket(NoisyFunction * f1d, NoisyFunctionValue &a, NoisyFunctionValue &b, NoisyFunctionValue &c)
   {
      using namespace std;

      if (f1d->getNDim() != 1) cout << "ERROR parabgoldMinimization(): The NoisyFunction is not 1D. Ndim=" << f1d->getNDim()  << endl;
      
      NoisyFunctionValue fooa(1), foob(1), fooc(1);

      //cout << "a=" << a.getX(0) << endl;
      b.setX(a.getX(0)+0.25);
      //cout << "b=" << b.getX(0) << endl;
      double newf, dnewf;
      f1d->f(b.getX(),newf,dnewf);
      b.setF(newf,dnewf);
      while(b==a)
      {
         //cout << "WHILE 1    fa=" << a.getF() << "   fb=" << b.getF() << endl;
         //cout << "dfa=" << a.getDf() << "   dfb=" << b.getDf() << endl;
         b.setX(a.getX(0)+1.5*(b.getX(0)-a.getX(0)));
         //cout << "b=" << b.getX(0) << endl;
         f1d->f(b.getX(),newf,dnewf);
         b.setF(newf,dnewf);
      }
      
      if (b<a)
      {
         fooa=a;
         foob=b;
         //cout << "COND b<a" << endl;
         c.setX(foob.getX(0)+(foob.getX(0)-fooa.getX(0)));
         //cout << "c=" << c.getX(0) << endl;
         f1d->f(c.getX(),newf,dnewf);
         c.setF(newf,dnewf);
         while (!(c>b)){
            if (c<foob)
            {
               a=b;
               b=c;
            }
            //cout << "WHILE 2" << endl;
            c.setX(foob.getX(0)+1.5*(c.getX(0)-foob.getX(0)));
            //cout << "c=" << c.getX(0) << endl;
            f1d->f(c.getX(),newf,dnewf);
            c.setF(newf,dnewf);
         }
      } else
      {
         fooa=a;
         foob=b;
         //cout << "COND b>a" << endl;
         c.setX(fooa.getX(0)-(foob.getX(0)-fooa.getX(0)));
         //cout << "c=" << c.getX(0) << endl;
         f1d->f(c.getX(),newf,dnewf);
         c.setF(newf,dnewf);
         while (!(c>a))
         {
            if (c<fooa)
            {
               b=a;
               a=c;
            }
            //cout << "WHILE 2" << endl;
            c.setX(fooa.getX(0)+1.5*(c.getX(0)-fooa.getX(0)));
            //cout << "c=" << c.getX(0) << endl;
            f1d->f(c.getX(),newf,dnewf);
            c.setF(newf,dnewf);
         }
         fooa=a; foob=b; fooc=c;
         a=fooc;
         b=fooa;
         c=foob;
      }

   }


   void parabgoldMinimization(NoisyFunction * f1d, const double &eps, NoisyFunctionValue &a, NoisyFunctionValue &b, NoisyFunctionValue &c)
   {
      using namespace std;

      if (f1d->getNDim() != 1) cout << "ERROR parabgoldMinimization(): The NoisyFunction is not 1D. Ndim=" << f1d->getNDim()  << endl;

      int count=0, lh=0, rh=0, cx=0;
      bool force_lh=false, force_rh=false;

      if (b>a) cout << "ERROR parabgoldMinimization(): b>a" << endl;
      if (b>c) cout << "ERROR parabgoldMinimization(): b>c" << endl;
      
      double present_eps=(a>c)?(a.getF()-b.getF()-a.getDf()-b.getDf()):(c.getF()-b.getF()-c.getDf()-b.getDf());
      
      NoisyFunctionValue x(1);
      double delta;
      double newf, dnewf;

      //cout << endl << endl << "a=" << a.getX(0) << "     b=" << b.getX(0) << "    c=" << c.getX(0) << endl;
      //cout << "fa=" << a.getF()  << "    fb=" << b.getF()  << "   fc=" << c.getF() << endl << endl;
      
      while (present_eps>eps)
      {
         //update counter
         count++;
         //check the conditions for a,b and c
         if (b>a) cout << "ERROR parabgoldMinimization(): b>a, count=" << count << endl;
         if (b>c) cout << "ERROR parabgoldMinimization(): b>c, count=" << count << endl;
         
         //find the new x using a parabolic interpolation
         x.setX( b.getX(0)-0.5*((b.getX(0)-a.getX(0))*(b.getX(0)-a.getX(0))*(b.getF()-c.getF())-(b.getX(0)-c.getX(0))*(b.getX(0)-c.getX(0))*(b.getF()-a.getF()))/((b.getX(0)-a.getX(0))*(b.getF()-c.getF())-(b.getX(0)-c.getX(0))*(b.getF()-a.getF())) );
         
         //if the range has been reduced on the same side of b too many times, 
         //or last time the parabolic prediction has led to a point equivalent to b, 
         //use the golden bi-section rule on the other side of b
         if (((x.getX(0)<b.getX(0))&&(lh>1))||((x.getX(0)>b.getX(0))&&(rh>1))||(cx>0))
         {
            if (b.getX(0)-a.getX(0)>c.getX(0)-b.getX(0))
            {
               delta=(b.getX(0)-a.getX(0))*0.38197*0.5;
               x.setX(b.getX(0)-delta);
            }
            else
            {
               delta=(c.getX(0)-b.getX(0))*0.38197*0.5;
               x.setX(b.getX(0)+delta);
            }
         }
         
         //if it has been asked to propose a point on the left side of b
         if (force_lh){
            delta=(b.getX(0)-a.getX(0))*0.38197*0.5;
            x.setX(b.getX(0)-delta);
            //cout << "force SX" << endl;
         }
         //if it has been asked to propose a point on the right side of b
         if (force_rh){
            delta=(c.getX(0)-b.getX(0))*0.38197*0.5;
            x.setX(b.getX(0)+delta);
            //cout << "force DX" << endl;
         }
         
         //compute the value of f in the point x
         f1d->f(x.getX(),newf,dnewf);
         x.setF(newf,dnewf);
         
         //cout << "x=" << x.getX(0) << "     b=" << b.getX(0) << endl;
         force_rh=false; force_lh=false;
         if (x.getX(0)<b.getX(0))
         {
            lh++; rh=0;
            if (x>b)
            {
               cx=0;
               //cout << "a=x" << endl;
               a=x;
            } else if (x<b)
            {
               cx=0;
               //cout << "c=b" << endl;
               c=b;
               //cout << "b=x" << endl;
               b=x;
            } else
            {
               if (cx>0) {force_rh=true; force_lh=false;}
               cx++;
               if (x.getF()<b.getF())
               {
                  //cout << "b=x" << endl;
                  b=x;
               }
            }
            
         } else
         {
            lh=0; rh++;
            if (x>b)
            {
               cx=0;
               //cout << "c=x" << endl;
               c=x;
            } else if (x<b)
            {
               cx=0;
               //cout << "a=b" << endl;
               a=b;
               //cout << "b=x" << endl;
               b=x;
            } else
            {
               if (cx>0) {force_lh=true; force_rh=false;}
               cx++;
               if (x.getF()<b.getF())
               {
                  //cout << "b=x" << endl;
                  b=x;
               }
            }
         }
         //update the actual eps value
         present_eps=(a>c)?(a.getF()-b.getF()-a.getDf()-b.getDf()):(c.getF()-b.getF()-c.getDf()-b.getDf());
         //cout << "a=" << a.getX(0) << "     b=" << b.getX(0) << "    c=" << c.getX(0) << endl;
         //cout << "fa=" << a.getF()  << "    fb=" << b.getF()  << "   fc=" << c.getF() << endl << endl;
         if (cx>2) break;
      }
      //cout << "ParabGold Terminated. count=" << count << endl;
   }



}



