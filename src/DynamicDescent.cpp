#include "nfm/DynamicDescent.hpp"

#include "nfm/LogManager.hpp"

#include <cmath>
#include <numeric>

namespace nfm
{

DynamicDescent::DynamicDescent(NoisyFunctionWithGradient * targetfun, const int max_n_const_values,
                               const double stepSize, const double alpha):
        NFM(targetfun, max_n_const_values), _stepSize(stepSize), _alpha(alpha)
{}

// --- Minimization

void DynamicDescent::_findMin()
{
    LogManager::logString("\nBegin DynamicDescent::findMin() procedure\n");

    // vector to hold the gradient and (possibly) error
    std::vector<NoisyValue> grad(static_cast<size_t>(_ndim));

    // holds previous old update (init to 0)
    std::vector<double> dx(grad.size()); // stores momentum updates

    //begin the minimization loop
    int cont = 0;
    while (true) {
        ++cont;
        LogManager::logString("\n\nDynamicDescent::findMin() Step " + std::to_string(cont) + "\n");

        // compute the gradient and current target
        _last.f = _gradfun->fgrad(_last.x, grad);
        this->_storeLastValue();
        this->_writeGradientToLog(grad);
        if (this->_shouldStop(&grad)) { break; }

        // find the next position
        this->_findNextX(grad, dx);
    }

    LogManager::logNoisyIOPair(_last, LogLevel::NORMAL, "Final position and target value");
    LogManager::logString("\nEnd DynamicDescent::findMin() procedure\n");
}

// --- Internal methods

void DynamicDescent::_findNextX(const std::vector<NoisyValue> &grad, std::vector<double> &dx)
{
    // update momenta and _last
    for (int i = 0; i < _ndim; ++i) {
        dx[i] = _alpha*dx[i] - _stepSize*grad[i].value;
        _last.x[i] += dx[i];
    }
    // report to the log
    this->_writeXUpdateToLog(dx);
}
} // namespace nfm