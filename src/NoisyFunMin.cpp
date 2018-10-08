#include "NoisyFunMin.hpp"
#include "LogNFM.hpp"

#include <cmath>
#include <iostream>
#include <sstream>

// --- Protected methods

bool NFM::_isNotConverged(){
    using namespace std;

    NoisyFunctionValue * v = new NoisyFunctionValue(_x->getNDim());
    *v = *_x;
    _old_values.push_front(v);

    if (_old_values.size() > _max_n_const_values) {
        delete _old_values.back();
        _old_values.pop_back();
    }

    if (_old_values.size() == _max_n_const_values) {
        for (list<NoisyFunctionValue *>::iterator it = _old_values.begin(); it != _old_values.end(); ++it){
            if (it != _old_values.begin()){
                if (! (**it == **_old_values.begin()) ){
                    return true;
                }
            }
        }
        if (_old_values.size() < _max_n_const_values) {
            return true;
        }

        NFMLogManager * log_manager = new NFMLogManager();
        log_manager->writeOnLog("\nCost function has stabilised, interrupting minimisation procedure.\n");
        delete log_manager;

        return false;
    }

    return true;
}


bool NFM::_meaningfulGradient(const double * grad, const double * graderr)
{
    for (int i=0; i<_ndim; ++i)
        {
            if (std::abs(grad[i])>graderr[i]) return true;
        }
    return false;
}


// --- Loggers

void NFM::_writeCurrentXInLog()
{
    NFMLogManager log_manager = NFMLogManager();
    log_manager.writeNoisyValueInLog(_x, "Current position and target value");
}

void NFM::_writeGradientInLog(const double * grad, const double * dgrad)
{
    NFMLogManager log_manager = NFMLogManager();
    log_manager.writeVectorInLog(grad, dgrad, _ndim, "Raw gradient", "g");
}

void NFM::_writeXUpdateInLog(const double * xu)
{
    NFMLogManager log_manager = NFMLogManager();
    log_manager.writeVectorInLog(xu, NULL, _ndim, "Position update", "u");
}

void NFM::_reportMeaninglessGradientInLog()
{
    using namespace std;

    NFMLogManager log_manager = NFMLogManager();

    stringstream s;
    s << endl << "gradient seems to be meaningless, i.e. its error is too large" << endl;
    s << flush;
    log_manager.writeOnLog(s.str());
}

void NFM::_writeOldValuesInLog()
{
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

    for (NoisyFunctionValue * v : _old_values) {
        delete v;
    }
        _old_values.clear();
}
