#ifndef ORNATE_MATH_VECTOR_H
#define ORNATE_MATH_VECTOR_H

#include "math_utils.h"

using std::vector;

namespace ornate {
template <typename T>
inline double l2_norm(const vector<T>& a, const vector<T>& b) {
    double ret = 0.0;
    size_t l = a.size();
    for (size_t i = 0; i < l; ++i) {
        ret += std::pow(a[i] - b[i], 2);
    }
    return ret;
}

/**
 * vs_*_inplace means vector op scalar, calc in place
 */
template <typename T, typename T1>
inline void vs_multiply_inplace(vector<T>& a, const T1 x) {
    size_t l = a.size();
    for (size_t i = 0; i < l; ++i) {
        a[i] *= x;
    }
}

template <typename T, typename T1>
inline void vs_add_inplace(vector<T>& a, const T1 x) {
    size_t l = a.size();
    for (size_t i = 0; i < l; ++i) {
        a[i] += x;
    }
}

template <typename T, typename T1>
inline void vs_minus_inplace(vector<T>& a, const T1 x) {
    size_t l = a.size();
    for (size_t i = 0; i < l; ++i) {
        a[i] -= x;
    }
}

template <typename T, typename T1>
inline void vs_divide_inplace(vector<T>& a, const T1 x) {
    size_t l = a.size();
    for (size_t i = 0; i < l; ++i) {
        a[i] /= x;
    }
}

/**
 * vs_* means vector op scalar, calc produce new vector
 */
template <typename T, typename T1, typename TRet = T>
inline vector<TRet> vs_multiply(const vector<T>& a, const T1 x) {
    vector<TRet> ret;
    size_t l = a.size();
    ret.resize(l);
    for (size_t i = 0; i < l; ++i) {
        ret[i] = a[i] * x;
    }
    return ret;
}
template <typename T, typename T1, typename TRet = T>
inline vector<TRet> vs_add(const vector<T>& a, const T1 x) {
    vector<TRet> ret;
    size_t l = a.size();
    ret.resize(l);
    for (size_t i = 0; i < l; ++i) {
        ret[i] = a[i] + x;
    }
    return ret;
}
template <typename T, typename T1, typename TRet = T>
inline vector<TRet> vs_minus(const vector<T>& a, const T1 x) {
    vector<TRet> ret;
    size_t l = a.size();
    ret.resize(l);
    for (size_t i = 0; i < l; ++i) {
        ret[i] = a[i] - x;
    }
    return ret;
}
template <typename T, typename T1, typename TRet = T>
inline vector<TRet> vs_divide(const vector<T>& a, const T1 x) {
    vector<TRet> ret;
    size_t l = a.size();
    ret.resize(l);
    for (size_t i = 0; i < l; ++i) {
        ret[i] = a[i] / x;
    }
    return ret;
}

/**
 * vv means vector vector element wise
 */
template <typename T, typename T1, typename TRet = T>
inline vector<TRet> vv_multiply(const vector<T>& a, const vector<T1>& b) {
    vector<TRet> ret;
    size_t l = a.size();
    ret.resize(l);
    for (size_t i = 0; i < l; ++i) {
        ret[i] = a[i] * b[i];
    }
    return ret;
}
template <typename T, typename T1, typename TRet = T>
inline vector<TRet> vv_add(const vector<T>& a, const vector<T1>& b) {
    vector<TRet> ret;
    size_t l = a.size();
    ret.resize(l);
    for (size_t i = 0; i < l; ++i) {
        ret[i] = a[i] + b[i];
    }
    return ret;
}
template <typename T, typename T1, typename TRet = T>
inline vector<TRet> vv_minus(const vector<T>& a, const vector<T1>& b) {
    vector<TRet> ret;
    size_t l = a.size();
    ret.resize(l);
    for (size_t i = 0; i < l; ++i) {
        ret[i] = a[i] - b[i];
    }
    return ret;
}
template <typename T, typename T1, typename TRet = T>
inline vector<TRet> vv_divide(const vector<T>& a, const vector<T1>& b) {
    vector<TRet> ret;
    size_t l = a.size();
    ret.resize(l);
    for (size_t i = 0; i < l; ++i) {
        ret[i] = a[i] / b[i];
    }
    return ret;
}
}  // namespace ornate

#endif
