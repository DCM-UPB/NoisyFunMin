#include "nfm/NoisyFunMin.hpp"
#include "nfm/LogManager.hpp"

#include <cmath>
#include <iostream>
#include <sstream>

namespace nfm
{

// --- Protected methods

void NFM::_clearOldValues()
{
    for (NoisyValue * v : _old_values) {
        delete v;
    }
    _old_values.clear();
}


void NFM::_storeOldValue()
{
    if (_max_n_const_values > 0) {
        auto * v = new NoisyValue(_x->getNDim());
        *v = *_x;
        _old_values.push_front(v);

        if (_old_values.size() > _max_n_const_values) {
            delete _old_values.back();
            _old_values.pop_back();
        }
    }
}


bool NFM::_isConverged()
{
    using namespace std;

    if (_max_n_const_values < 1) { return false; }

    if (_old_values.size() == _max_n_const_values) {
        for (auto it = _old_values.begin(); it != _old_values.end(); ++it) {
            if (it != _old_values.begin()) {
                if (!(**it == **_old_values.begin())) {
                    return false;
                }
            }
        }

        LogManager log_manager;
        log_manager.writeOnLog("\nCost function has stabilised, interrupting minimization procedure.\n");

        return true;
    }

    return false;
}


bool NFM::_meaningfulGradient(const double * grad, const double * graderr)
{
    if (_useGradientError) {
        for (int i = 0; i < _ndim; ++i) {
            if (fabs(grad[i]) > graderr[i]) { return true; }
        }
    }
    else {
        return true;
    }

    LogManager log_manager;
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
    LogManager log_manager;
    if (log_manager.isVerbose()) {
        log_manager.writeNoisyValueInLog(_x, 2, "Current position and target value", "f", true, "x");
    }
    else { log_manager.writeNoisyValueInLog(_x, 1, "Current target value", "f", false); }
}

void NFM::_writeGradientInLog(const double * grad, const double * dgrad)
{
    LogManager log_manager;
    log_manager.writeVectorInLog(grad, _useGradientError ? dgrad : nullptr, _ndim, 2, "Raw gradient", "g");
}

void NFM::_writeXUpdateInLog(const double * xu)
{
    LogManager log_manager;
    log_manager.writeVectorInLog(xu, nullptr, _ndim, 2, "Position update", "u");
}

void NFM::_writeOldValuesInLog()
{
    using namespace std;

    LogManager log_manager;

    stringstream s;
    s << endl << "last values:    ";
    for (auto &_old_value : _old_values) {
        s << _old_value->getF() << " +- " << _old_value->getDf() << "    ";
    }
    s << endl;
    s << "equal to first element? ";
    for (auto it = _old_values.begin(); it != _old_values.end(); ++it) {
        if (it != _old_values.begin()) {
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
    _x = new NoisyValue(_ndim);
    for (int i = 0; i < _ndim; ++i) { _x->setX(i, 0.); }
    //gradient of the target function
    _gradtargetfun = nullptr;
    _flaggradtargetfun = false;
    //optimization's domain
    //_indomain = 0;
    //_flagindomain = false;
    //eps
    _epstargetfun = 0.;
    _epsx = 0.;
}


NFM::~NFM()
{
    //deallocate x
    delete _x;

    _clearOldValues();
}
} // namespace nfm