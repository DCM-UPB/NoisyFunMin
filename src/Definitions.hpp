#ifndef DEFINITIONS
#define DEFINITIONS


namespace nfm
{
   //1D Function
   typedef void (*Fun1D)(const double &, double &, double &);
   //                       ^input          ^output   ^error on the output
   
   //                       ^input          ^output   ^error on the output
   typedef void (*TargetFun)(const int &, const double *, double &, double &);
   //                           ^dim of input  ^input        ^output   ^error on the output
   
   // Gradient of the target Function to minimize
   typedef void (*GradTargetFun)(const int &, const double *, double *, double *);
   //                              ^dim of input/output  ^input    ^output  ^error on the output
   
   // Domain of optimization: this function returns true if the point is in the domain, false otherwise
   typedef bool (*DomainFun)(const int &, const double *);
   //                        ^dim of input    ^input
}


#endif

