## perf
Handwritten code can still be faster than the Rcpp sugar code by a small factor

# Implementation
Rcpp sugar employ a technique called the Curiously Recurring Template Pattern (CRTP)

CRTP allows the base class to access methods of the derived class

template <int RTYPE, bool na, typename VECTOR> class VectorBase