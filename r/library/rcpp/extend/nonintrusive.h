#include <RcppCommon.h>

// third party library that declares class Bar
#include <foobar.h>

// declaring the specialization
namespace Rcpp {
template <>
SEXP wrap(const Bar&);
template <>
Bar as(SEXP) throw(not_compatible);

namespace traits {
template <typename T>
SEXP wrap(const Bling<T>&);

template <typename T>
class Exporter {
public:
    Exporter(SEXP x) : t(x) {}
    inline T get() { return t; }

private:
    T t;
};

template <typename T>
class Exporter<Bling<T> >;
}  // namespace traits
}  // namespace Rcpp

// this must appear after the specialization, else the specialization will not be seen by Rcpp types
#include <Rcpp.h>
