#ifndef ORNATE_MATH_MATRIX_H
#define ORNATE_MATH_MATRIX_H

#include "math_utils.h"

using std::vector;

namespace ornate {
template <typename T>
inline vector<vector<T>> mm_multiply(const vector<vector<T>>& A, const vector<vector<T>>& vs) {
    size_t k = A.size(), d = A[0].size(), l = vs[0].size();
    vector<vector<T>> ret;
    ret.resize(k);
    for (size_t i = 0; i < k; ++i) {
        ret[i].resize(l);
        for (size_t j = 0; j < l; ++j) {
            T sum = 0;
            for (size_t m = 0; m < d; ++m) {
                sum += A[i][m] * vs[m][j];
            }
            ret[i][j] = sum;
        }
    }
    return ret;
}

template <typename T>
inline vector<vector<T>> transpose(const vector<vector<T>>& A) {
    size_t k = A.size(), d = A[0].size();
    vector<vector<T>> ret;
    ret.resize(d);
    for (size_t i = 0; i < d; ++i) {
        ret[i].resize(k);
        for (size_t j = 0; j < k; ++j) {
            ret[i][j] = A[j][i];
        }
    }
    return ret;
}
}  // namespace ornate

#endif
