#ifndef ORNATE_MATH_VECTOR_H
#define ORNATE_MATH_VECTOR_H

#include "math_helper.h"
#include "math_utils.h"

using std::vector;

namespace ornate {

namespace detail {
template <typename T, typename T1, template <typename, typename, typename> class TFunctor>
inline void __vs_op_inplace(vector<T>& a, const T1 x, const TFunctor<T, T1, T>& functor) {
    size_t l = a.size();
    for (size_t i = 0; i < l; ++i) {
        a[i] = functor(a[i], x);
    }
}

template <typename T, typename T1, typename TRet = T, template <typename, typename, typename> class TFunctor>
inline vector<TRet> __vs_op(const vector<T>& a, const T1 x, const TFunctor<T, T1, TRet>& functor) {
    vector<TRet> ret;
    size_t l = a.size();
    ret.resize(l);
    for (size_t i = 0; i < l; ++i) {
        ret[i] = functor(a[i], x);
    }
    return ret;
}

template <typename T, typename T1, typename TRet = T, template <typename, typename, typename> class TFunctor>
inline vector<TRet> __vv_op(const vector<T>& a, const vector<T1>& b, const TFunctor<T, T1, TRet>& functor) {
    vector<TRet> ret;
    size_t l = a.size();
    ret.resize(l);
    for (size_t i = 0; i < l; ++i) {
        ret[i] = functor(a[i], b[i]);
    }
    return ret;
}
}  // namespace detail
template <typename T>
inline double l2_norm(const vector<T>& a, const vector<T>& b) {
    double ret = 0.0;
    size_t l = a.size();
    for (size_t i = 0; i < l; ++i) {
        ret += std::pow(a[i] - b[i], 2);
    }
    return ret;
}

template <typename T>
std::vector<T> linspace(T a, T b, size_t N) {
    T h = (b - a) / static_cast<T>(N - 1);
    std::vector<T> xs(N);
    typename std::vector<T>::iterator x;
    T val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h) *x = val;
    return xs;
}

/**
 * vs_*_inplace means vector op scalar, calc in place
 */
template <typename T, typename T1>
inline void vs_multiply_inplace(vector<T>& a, const T1 x) {
    detail::__vs_op_inplace(a, x, multiply_dt<T, T1, T>());
}
template <typename T, typename T1>
inline void vs_add_inplace(vector<T>& a, const T1 x) {
    detail::__vs_op_inplace(a, x, plus_dt<T, T1, T>());
}
template <typename T, typename T1>
inline void vs_minus_inplace(vector<T>& a, const T1 x) {
    detail::__vs_op_inplace(a, x, minus_dt<T, T1, T>());
}
template <typename T, typename T1>
inline void vs_divide_inplace(vector<T>& a, const T1 x) {
    detail::__vs_op_inplace(a, x, divide_dt<T, T1, T>());
}

/**
 * vs_* means vector op scalar, calc produce new vector
 */
template <typename T, typename T1, typename TRet = T>
inline vector<TRet> vs_multiply(const vector<T>& a, const T1 x) {
    return detail::__vs_op(a, x, multiply_dt<T, T1, TRet>());
}
template <typename T, typename T1, typename TRet = T>
inline vector<TRet> vs_add(const vector<T>& a, const T1 x) {
    return detail::__vs_op(a, x, plus_dt<T, T1, TRet>());
}
template <typename T, typename T1, typename TRet = T>
inline vector<TRet> vs_minus(const vector<T>& a, const T1 x) {
    return detail::__vs_op(a, x, minus_dt<T, T1, TRet>());
}
template <typename T, typename T1, typename TRet = T>
inline vector<TRet> vs_divide(const vector<T>& a, const T1 x) {
    return detail::__vs_op(a, x, divide_dt<T, T1, TRet>());
}

/**
 * vv means vector vector element wise
 */
template <typename T, typename T1, typename TRet = T>
inline vector<TRet> vv_multiply(const vector<T>& a, const vector<T1>& b) {
    return detail::__vv_op(a, b, multiply_dt<T, T1, TRet>());
}
template <typename T, typename T1, typename TRet = T>
inline vector<TRet> vv_add(const vector<T>& a, const vector<T1>& b) {
    return detail::__vv_op(a, b, plus_dt<T, T1, TRet>());
}
template <typename T, typename T1, typename TRet = T>
inline vector<TRet> vv_minus(const vector<T>& a, const vector<T1>& b) {
    return detail::__vv_op(a, b, minus_dt<T, T1, TRet>());
}
template <typename T, typename T1, typename TRet = T>
inline vector<TRet> vv_divide(const vector<T>& a, const vector<T1>& b) {
    return detail::__vv_op(a, b, divide_dt<T, T1, TRet>());
}
}  // namespace ornate

#endif
