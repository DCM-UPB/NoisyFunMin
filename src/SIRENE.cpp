#include "nfm/SIRENE.hpp"

#include "nfm/LogManager.hpp"

#include <numeric>
#include <cmath>

#include <iostream>
namespace nfm
{

// --- Minimization

void SIRENE::_findMin()
{
    LogManager::logString("\nBegin SIRENE::findMin() procedure\n");

    // helper variables
    const std::vector<double> &F = _grad.val; // convenience reference (the force)
    const std::vector<double> &Ferr = _grad.err; // the error of force
    std::vector<double> v(F.size()); // velocity vector
    std::vector<double> m(F.size()); // gradient running average (only for beta>0)

    double dt = _dt0; // current time-step
    double alpha = _alpha0; // current mixing factor
    int Npos = 0; // number of steps since "m.v" was negative

    //begin the minimization loop
    int iter = 0;
    while (true) {
        ++iter;
        if (LogManager::isLoggingOn()) { // else skip string construction
            LogManager::logString("\nSIRENE::findMin() Step " + std::to_string(iter) + "\n");
        }
        std::cout << std::endl << "SIRENE step " << iter << ":" << std::endl;

        // compute the gradient and current target
        this->_updateTarget();
        if (this->_shouldStop()) { break; } // we are done

        // update vector lengths
        const double vnorm = sqrt(std::inner_product(v.begin(), v.end(), v.begin(), 0.));
        const double fnorm = sqrt(std::inner_product(F.begin(), F.end(), F.begin(), 0.));

        // calculate scalar product P = F.v or P = m.v (if beta > 0)
        /*        double P;
                if (_beta > 0) {
                    for (int i = 0; i < _ndim; ++i) {
                        m[i] = _beta*m[i] + (1. - _beta)*F[i];
                    }
                    const double mnorm = sqrt(std::inner_product(m.begin(), m.end(), m.begin(), 0.));
                    P = std::inner_product(m.begin(), m.end(), v.begin(), 0.);
                    P = (P != 0.) ? P/(mnorm*vnorm) : 0.; // P would be 0 if mnorm or vnorm were 0
                }
                else {
                    P = std::inner_product(F.begin(), F.end(), v.begin(), 0.);
                    P = (P != 0.) ? P/(fnorm*vnorm) : 0.; // P would be 0 if fnorm or vnorm were 0
                }
        */

        NoisyValue P{}; // in this algorithm P is a NoisyValue
        P.val = std::inner_product(F.begin(), F.end(), v.begin(), 0.);
        // make an estimation of the error based on error propagation
        for (int i = 0; i < _ndim; ++i) {
            P.err += pow(Ferr[i]*v[i], 2);
        }
        P.err = sqrt(P.err);

        /*if (P.val != 0.) {
            P /= (fnorm*vnorm);
        }*/
        std::cout << "vnorm " << vnorm << ", fnorm " << fnorm << ", P" << P << std::endl;

        // modify velocity
        for (int i = 0; i < _ndim; ++i) {
            v[i] = (1. - alpha)*v[i] + alpha*vnorm*F[i]/fnorm;
        }

        // check P
        if (P > 0.) { // we are going downhill
            std::cout << "P > 0" << std::endl;
            if (++Npos > _Nmin) { // then increase dt
                std::cout << "Npos(" << Npos << ") > Nmin(" << _Nmin << ")" << std::endl;
                dt = std::min(dt*_finc, _dtmax);
                alpha *= _falpha;
                std::cout << "dt " << dt << std::endl;
            }
        }
        else if (P < 0.) { // we are going uphill
            std::cout << "P < 0" << std::endl;
            Npos = 0;
            dt *= _fdec;
            alpha = _alpha0;
            if (_flag_soft) { // freeze softly
                for (int i=0; i<_ndim; ++i) {
                    const double p_err = fabs(NoisyValue::getSigmaLevel()*Ferr[i]*v[i]); // err on Fi*vi
                    if (F[i]*v[i] < -p_err) {
                        const double vold = v[i];
                        //v[i] = 0.;
                        v[i] = (F[i] != 0.) ? v[i] = p_err/F[i] : 0.;
                        std::cout << "flip v" << i << ": " << vold << " -> " << v[i] << std::endl;
                    }
                }
            }
            else { // freeze completely (default)
                std::fill(v.begin(), v.end(), 0.);
            }
            std::cout << "dt " << dt << std::endl;
        }
        else {
            //++Npos;
            std::cout << "P == 0" << std::endl;
            std::cout << "dt " << dt << std::endl;
        }

        // make MD step to find the next position/velocity
        this->_doMDStep(v, dt);
    }

    LogManager::logString("\nEnd SIRENE::findMin() procedure\n");
}
} // namespace nfm