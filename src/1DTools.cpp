#include "nfm/1DTools.hpp"

#include "nfm/LogNFM.hpp"

#include <iostream>
#include <sstream>
#include <stdexcept>


namespace nfm
{
    const double HUGE = 1.e+6;
    const int MAX_NUM_EVAL_FOR_BRACKET = 100;


    void findBracket(NoisyFunction * f1d, NoisyFunctionValue &a, NoisyFunctionValue &b, NoisyFunctionValue &c)
    {
        using namespace std;

        if (f1d->getNDim() != 1) {
            cout << "ERROR parabgoldMinimization(): The NoisyFunction is not 1D. Ndim=" << f1d->getNDim()  << endl;
        }

        double newf, dnewf;
        int count_newf = 0;

        //cout << "a=" << a.getX(0) << endl;
        b.setX(a.getX(0)+0.25);
        //cout << "b=" << b.getX(0) << endl;
        if (++count_newf > MAX_NUM_EVAL_FOR_BRACKET) { _abortFindBracket(); }

        f1d->f(b.getX(),newf,dnewf);
        b.setF(newf,dnewf);
        _writeabcInLog("findBracket init", a, b, b);

        while(b==a){   // if fb==fa, increase b until when the two values are different
            //cout << "WHILE 1    fa=" << a.getF() << "   fb=" << b.getF() << endl;
            //cout << "dfa=" << a.getDf() << "   dfb=" << b.getDf() << endl;
            b.setX(a.getX(0)+1.5*(b.getX(0)-a.getX(0)));
            //cout << "b=" << b.getX(0) << endl;
            if (++count_newf > MAX_NUM_EVAL_FOR_BRACKET) { _abortFindBracket(); }

            f1d->f(b.getX(),newf,dnewf);
            b.setF(newf,dnewf);
            _writeabcInLog("findBracket init", a, b, b);

            if (b.getX(0) > HUGE) {
                break;
            }
        }

        if (b==a){  // could not find a b such that fb!=fa, now try looking in the other direction
            b.setX(a.getX(0));
            b.setF(a.getF(), a.getDf());
            a.setX(b.getX(0)-0.25);
            f1d->f(a.getX(),newf,dnewf);
            if (++count_newf > MAX_NUM_EVAL_FOR_BRACKET) { _abortFindBracket(); }
            a.setF(newf,dnewf);
            _writeabcInLog("findBracket init", a, a, b);

            while(b==a){   // if fb==fa, increase b until when the two values are different
                //cout << "WHILE 1    fa=" << a.getF() << "   fb=" << b.getF() << endl;
                //cout << "dfa=" << a.getDf() << "   dfb=" << b.getDf() << endl;
                a.setX(b.getX(0)-1.5*(b.getX(0)-a.getX(0)));
                //cout << "b=" << b.getX(0) << endl;
                if (++count_newf > MAX_NUM_EVAL_FOR_BRACKET) { _abortFindBracket(); }

                f1d->f(a.getX(), newf, dnewf);
                a.setF(newf,dnewf);
                _writeabcInLog("findBracket init", a, a, b);

                if (a.getX(0) < -HUGE) {
                    throw std::runtime_error( "NoisyFunctionMin Error! Bracketing impossible. Cannot find a b such that fa!=fb" );
                }
            }
        }

        if (b<a)
            {
                //cout << "COND b<a" << endl;
                c.setX(b.getX(0)+(b.getX(0)-a.getX(0)));
                //cout << "c=" << c.getX(0) << endl;
                if (++count_newf > MAX_NUM_EVAL_FOR_BRACKET) { _abortFindBracket(); }

                f1d->f(c.getX(),newf,dnewf);
                c.setF(newf,dnewf);

                _writeabcInLog("findBracket (A)", a, b, c);

                while (!(c>b)){
                    //cout << "WHILE 2" << endl;
                    if (c<b)
                        {
                            a=b;
                            b=c;
                            c.setX(b.getX(0)+2.*(b.getX(0)-a.getX(0)));
                        } else {
                        c.setX(c.getX(0)+2.*(c.getX(0)-b.getX(0)));
                    }
                    //cout << "c=" << c.getX(0) << endl;
                    if (++count_newf > MAX_NUM_EVAL_FOR_BRACKET) { _abortFindBracket(); }
                    f1d->f(c.getX(), newf, dnewf);
                    c.setF(newf,dnewf);

                    _writeabcInLog("findBracket (B)", a, b, c);
                }
            } else
            {
                c=b;
                b=a;
                //cout << "COND b>a" << endl;
                a.setX(b.getX(0)-(c.getX(0)-b.getX(0)));
                //cout << "c=" << c.getX(0) << endl;
                if (++count_newf > MAX_NUM_EVAL_FOR_BRACKET) { _abortFindBracket(); }
                f1d->f(a.getX(), newf, dnewf);
                a.setF(newf, dnewf);

                _writeabcInLog("findBracket (C)", a, b, c);

                while (!(a>b))
                    {
                        //cout << "WHILE 2" << endl;
                        if (a<b)
                            {
                                c=b;
                                b=a;
                                a.setX(b.getX(0)-2.*(c.getX(0)-b.getX(0)));
                            } else {
                            a.setX(a.getX(0)-2.*(b.getX(0)-a.getX(0)));
                        }
                        //cout << "c=" << c.getX(0) << endl;
                        if (++count_newf > MAX_NUM_EVAL_FOR_BRACKET) { _abortFindBracket(); }
                        f1d->f(a.getX(), newf, dnewf);
                        a.setF(newf, dnewf);

                        _writeabcInLog("findBracket (D)", a, b, c);
                    }
            }

        _writeabcInLog("findBracket", a, b, c);
    }


    void parabgoldMinimization(NoisyFunction * f1d, const double &eps, NoisyFunctionValue &a, NoisyFunctionValue &b, NoisyFunctionValue &c)
    {
        using namespace std;

        if (f1d->getNDim() != 1) {
            cout << "ERROR parabgoldMinimization(): The NoisyFunction is not 1D. Ndim=" << f1d->getNDim()  << endl;
        }

        int count=0, lh=0, rh=0, cx=0;
        bool force_lh=false, force_rh=false;

        if (b>a) { cout << "ERROR parabgoldMinimization(): b>a" << endl; }
        if (b>c) { cout << "ERROR parabgoldMinimization(): b>c" << endl; }

        double present_eps=(a>c)?(a.getF()-b.getF()-a.getDf()-b.getDf()):(c.getF()-b.getF()-c.getDf()-b.getDf());

        NoisyFunctionValue x(1);
        double delta;
        double newf, dnewf;

        //cout << endl << endl << "a=" << a.getX(0) << "     b=" << b.getX(0) << "    c=" << c.getX(0) << endl;
        //cout << "fa=" << a.getF()  << "    fb=" << b.getF()  << "   fc=" << c.getF() << endl << endl;

        _writeabcInLog("parabgoldMinimization", a, b, c);

        while (present_eps>eps)
            {
                //update counter
                count++;
                //check the conditions for a,b and c
                if (b>a) { cout << "ERROR parabgoldMinimization(): b>a, count=" << count << endl; }
                if (b>c) { cout << "ERROR parabgoldMinimization(): b>c, count=" << count << endl; }

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
                    }
                else
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
                if (cx>2) { break; }
            }
        //cout << "ParabGold Terminated. count=" << count << endl;

        _writeabcInLog("parabgoldMinimization", a, b, c);
    }



    void _writeabcInLog(const std::string &key, NoisyFunctionValue &a, NoisyFunctionValue &b, NoisyFunctionValue &c){
        using namespace std;

        NFMLogManager log_manager;

        stringstream s;
        s << key << ":    " <<
            a.getX(0) << " -> " << a.getF() << "    " <<
            b.getX(0) << " -> " << b.getF() << "    " <<
            c.getX(0) << " -> " << c.getF();
        s << flush;
        log_manager.writeOnLog(s.str(), 2);
    }



    void _abortFindBracket(){
        throw std::runtime_error( "NoisyFunctionMin Error! Bracketing is taking way too long..." );
    }



} // namespace nfm
