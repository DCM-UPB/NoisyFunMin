#include "nfm/DynamicDescent.hpp"

#include "nfm/LogManager.hpp"

#include <cmath>
#include <numeric>

namespace nfm
{

DynamicDescent::DynamicDescent(NoisyFunctionWithGradient * targetfun, const double stepSize, const int max_n_const_values):
        NFM(targetfun, max_n_const_values), _stepSize(stepSize)
{
    _inertia.assign(static_cast<size_t>(_ndim), 1.);
    _old_norm_dir.assign(static_cast<size_t>(_ndim), 0.);
}


// --- Log

void DynamicDescent::_writeInertiaToLog()
{
    LogManager::logVector(_inertia, LogLevel::VERBOSE, "Current inertia", "i");
}


// --- Minimization

void DynamicDescent::findMin()
{
    LogManager::logString("\nBegin DynamicDescent::findMin() procedure\n");

    // clear old values
    this->_clearOldValues();

    //arrays to hold the gradient and (possibly) error
    std::vector<NoisyValue> grad(static_cast<size_t>(_ndim));

    // init inertia to 1
    for (double &ini : _inertia) { ini = 1.; }

    //begin the minimization loop
    int cont = 0;
    while (true) {
        // compute the gradient and current target
        _last.f = _gradfun->fgrad(_last.x, grad);
        this->_storeOldValue();
        if (this->_shouldStop(&grad)) { break; }

        LogManager::logString("\n\nDynamicDescent::findMin() Step " + std::to_string(cont + 1) + "\n");
        this->_writeCurrentXToLog();
        this->_writeGradientToLog(grad);

        // find the next position
        this->findNextX(grad);

        cont++;
    }

    LogManager::logNoisyIOPair(_last, LogLevel::NORMAL, "Final position and target value");
    LogManager::logString("\nEnd DynamicDescent::findMin() procedure\n");
}

// --- Internal methods

void DynamicDescent::findNextX(const std::vector<NoisyValue> &grad)
{
    // compute the normalized gradection vector
    std::vector<double> norm_grad(static_cast<size_t>(_ndim));
    for (int i = 0; i< _ndim; ++i) { norm_grad[i] = grad[i].value; }
    const double gsum = sqrt(std::inner_product(norm_grad.begin(), norm_grad.end(), norm_grad.begin(), 0.0));
    if (gsum != 0.) {
        for (double &ngi : norm_grad) { ngi /= gsum; }
    }
    else {
        for (double &ngi : norm_grad) { ngi = 0.; }
    }

    // update the inertia
    for (int i = 0; i < _ndim; ++i) {
        const double iniUpdate = std::max(0., _old_norm_dir[i]*norm_grad[i]); // > 0 if the direction stays similar
        _inertia[i] *= 0.9 + 0.1*_ndim*iniUpdate; // exponential update of normalized inertia
    }
    // report it in the log
    this->_writeInertiaToLog();

    // update _last
    std::vector<double> dx(static_cast<size_t>(_ndim));
    for (int i = 0; i < _ndim; ++i) {
        dx[i] = -_stepSize*_inertia[i]*grad[i].value;
        _last.x[i] += dx[i];
    }
    this->_writeXUpdateToLog(dx);

    // store the grad for the next iteration
    _old_norm_dir = norm_grad;
}
} // namespace nfm