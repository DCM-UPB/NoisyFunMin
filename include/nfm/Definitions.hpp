#ifndef NFM_DEFINITIONS_HPP
#define NFM_DEFINITIONS_HPP


namespace nfm
{
    //1D Function
    using Fun1D = void (*)(const double &, double &, double &);
    //                       ^input          ^output   ^error on the output

    //                       ^input          ^output   ^error on the output
    using TargetFun = void (*)(const int &, const double *, double &, double &);
    //                           ^dim of input  ^input        ^output   ^error on the output

    // Gradient of the target Function to minimize
    using GradTargetFun = void (*)(const int &, const double *, double *, double *);
    //                              ^dim of input/output  ^input    ^output  ^error on the output

    // Domain of optimization: this function returns true if the point is in the domain, false otherwise
    using DomainFun = bool (*)(const int &, const double *);
    //                        ^dim of input    ^input
} // namespace nfm


#endif
