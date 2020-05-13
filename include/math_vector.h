#ifndef ORNATE_MATH_VECTOR_H
#define ORNATE_MATH_VECTOR_H

#include <algorithm>
#include <cstring>
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

template <typename T>
class __SortHelper {
public:
    const T* data{nullptr};
    int index;
    __SortHelper() = default;
    __SortHelper(const T& data_, int index_) : data{&data_}, index{index_} {}
    bool operator<(const __SortHelper<T>& a) const { return (*data < *a.data); }
};
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
 * get a slice from a vector. [from, to) can be < 0
 * behave much like in python, except that when to==0, slice to the end
 */
template <typename T>
std::vector<T> slice(const std::vector<T>& data, int from = 0, int to = 0) {
    if (data.empty()) return {};
    int total = static_cast<int>(data.size());
    int start_idx = from;
    int to_idx = to - 1;
    if (start_idx < 0) start_idx += total;
    if (to_idx < 0) to_idx += total;
    if (to_idx >= total) to_idx = total - 1;
    if (start_idx < 0 || start_idx >= total || to_idx < 0 || to_idx < start_idx) return {};
    std::vector<T> ret(static_cast<size_t>(to_idx - start_idx + 1));
    memcpy(&ret[0], &data[start_idx], sizeof(T) * ret.size());
    return ret;
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

template <typename T>
std::vector<int> get_sorted_index(const std::vector<T>& data, bool ascending_order = true) {
    int n = static_cast<int>(data.size());
    std::vector<detail::__SortHelper<T> > s(n);
    for (int i = 0; i < n; i++) {
        s[i].data = &data[i];
        s[i].index = i;
    }
    sort(s.begin(), s.end());
    std::vector<int> ret(n);
    for (int i = 0; i < n; i++) {
        ret[i] = s[i].index;
    }
    if (!ascending_order) {
        std::reverse(ret.begin(), ret.end());
    }
    return ret;
}

template <typename T>
std::vector<int> get_sorted_rank(const std::vector<T>& data, bool ascending_order = true) {
    int n = static_cast<int>(data.size());
    std::vector<detail::__SortHelper<T> > s(n);
    for (int i = 0; i < n; i++) {
        s[i].data = &data[i];
        s[i].index = i;
    }
    sort(s.begin(), s.end());
    std::vector<int> ret(n);
    for (int i = 0; i < n; i++) {
        ret[s[i].index] = ascending_order ? i : n - 1 - i;
    }
    return ret;
}

/**
 * filter out nan IN PLACE, move non-NAN to left
 * @param result_index, indices of non-NAN elements is saved if it is not nullptr
 * @return the number of non-NAN elements
 */
template <typename T>
uint32_t filter(INOUT T* n, uint32_t num_of_ele, OUT uint32_t* result_index = nullptr) {
    uint32_t c = 0;
    for (uint32_t i = 0; i < num_of_ele; i++) {
        if (isvalid(n[i])) {
            n[c] = n[i];
            if (result_index) result_index[c] = i;
            c++;
        }
    }
    std::fill(n + c, n + num_of_ele, NAN);
    return c;
}

/**
 * filter out nan IN PLACE. return the index of non-NAN elements
 */
template <typename T>
std::vector<uint32_t> filter(INOUT std::vector<T>& n) {
    std::vector<uint32_t> ret;
    if (n.size() == 0) return ret;
    uint32_t* index = new uint32_t[n.size()];
    uint32_t num = filter(&n[0], n.size(), index);
    ret.assign(index, index + num);
    n.resize(num);
    delete[] index;
    return ret;
}

/**
 * rank in place. will filter NAN, normalize to [0, 1], NAN will remain NAN
 */
template <typename T>
void rank(INOUT std::vector<T>& n) {
    std::vector<T> n2 = n;
    std::vector<uint32_t> index = filter(n2);
    std::vector<int> sorted = get_sorted_index(n2);
    std::vector<T> output(n.size(), NAN);
    int s = static_cast<int>(n2.size());
    if (s >= 2) {
        int start = 0;
        while (start < s) {
            int end_ = start + 1;
            auto curr_val = n2[sorted[start]];
            while (end_ < s){
                if(std::fabs(n2[sorted[end_]] - curr_val) > 1e-14) break;
                ++end_;
            }
            for (int jj = start; jj < end_; ++jj) {
                output[index[sorted[jj]]] = (start + end_ - 1) / (2.0 * (s - 1));
            }
            start = end_;
        }
    }
    n = output;
}

template <typename T>
float mean(const T* data, int32_t n) {
    T ret = 0;
    uint32_t count = 0;
    for (int32_t i = 0; i < n; i++) {
        if (isvalid(data[i])) {
            ret += data[i];
            count++;
        }
    }
    if (count > 0)
        return ret / count;
    else
        return NAN;
}

/**
 * get mean. will consider nan, [start_idx, end_idx), -1 to use all
 */
template <typename T>
float mean(IN const std::vector<T>& n, int32_t start_idx = -1, int32_t end_idx = -1) {
    if (start_idx < 0) start_idx = 0;
    if (end_idx < 0) end_idx = static_cast<int32_t>(n.size());
    return mean(n.data(), end_idx - start_idx);
}

/**
 * get sum. will consider nan.  [start_idx, end_idx), -1 to use all
 */
template <typename T>
float sum(IN const std::vector<T>& n, int32_t start_idx = -1, int32_t end_idx = -1) {
    if (start_idx < 0) start_idx = 0;
    if (end_idx < 0) end_idx = static_cast<int32_t>(n.size());
    float ret = 0;
    for (int32_t i = start_idx; i < end_idx; i++) {
        if (isvalid(n[i])) ret += n[i];
    }
    return ret;
}

/**
 * normalize. will consider nan
 */
template <typename T>
inline void normalize(INOUT std::vector<T>& n) {
    size_t s = n.size();
    float m = mean(n);
    for (size_t ii = 0; ii < s; ii++)
        if (isvalid(n[ii])) n[ii] -= m;
}

template <typename T>
inline void powerf(INOUT std::vector<T>& n, T exp) {
    size_t s = n.size();
    for (size_t ii = 0; ii < s; ii++)
        if (isvalid(n[ii])) {
            if (n[ii] < 0) n[ii] = -std::pow(-n[ii], exp);
            if (n[ii] > 0) n[ii] = std::pow(n[ii], exp);
        }
}

template <typename T>
inline void power(INOUT std::vector<T>& n, T exp, bool rank_first) {
    if (rank_first) {
        ornate::rank(n);
        normalize(n);
    }
    powerf(n, exp);
}
}  // namespace ornate

#endif
