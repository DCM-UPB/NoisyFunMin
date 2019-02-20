#include "nfm/NoisyFunMin.hpp"
#include "nfm/LogNFM.hpp"

#include <cmath>
#include <iostream>
#include <sstream>

// --- Protected methods

void NFM::_clearOldValues()
{
    for (NoisyFunctionValue * v : _old_values) {
        delete v;
    }
    _old_values.clear();
}

bool NFM::_isConverged(){
    using namespace std;

    if (_max_n_const_values < 1) return false;

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
                    return false;
                }
            }
        }

        NFMLogManager log_manager = NFMLogManager();
        log_manager.writeOnLog("\nCost function has stabilised, interrupting minimization procedure.\n");

        return true;
    }

    return false;
}


bool NFM::_meaningfulGradient(const double * grad, const double * graderr)
{
    if (_useGradientError) {
        for (int i=0; i<_ndim; ++i) {
            if (fabs(grad[i])>graderr[i]) return true;
        }
    }
    else {
        return true;
    }

    NFMLogManager log_manager = NFMLogManager();
    log_manager.writeOnLog("\nGradient seems to be meaningless, i.e. its error is too large.\n");
    return false;
}

bool NFM::_shouldStop(const double * grad, const double * graderr)
{
    return _isConverged() || !_meaningfulGradient(grad, graderr);
}


// --- Loggers

void NFM::_writeCurrentXInLog()
{
    NFMLogManager log_manager = NFMLogManager();
    if (log_manager.isVerbose()) log_manager.writeNoisyValueInLog(_x, 2, "Current position and target value", "f", true, "x");
    else log_manager.writeNoisyValueInLog(_x, 1, "Current target value", "f", false);
}

void NFM::_writeGradientInLog(const double * grad, const double * dgrad)
{
    NFMLogManager log_manager = NFMLogManager();
    log_manager.writeVectorInLog(grad, _useGradientError ? dgrad : NULL, _ndim, 2, "Raw gradient", "g");
}

void NFM::_writeXUpdateInLog(const double * xu)
{
    NFMLogManager log_manager = NFMLogManager();
    log_manager.writeVectorInLog(xu, NULL, _ndim, 2, "Position update", "u");
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

NFM::NFM(NoisyFunction * targetfun, const bool useGradientError, const size_t &max_n_const_values)
    : _useGradientError(useGradientError), _max_n_const_values(max_n_const_values)
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

    _clearOldValues();
}
