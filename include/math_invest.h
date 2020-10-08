#ifndef ORNATE_MATH_INVEST_H
#define ORNATE_MATH_INVEST_H

#include <deque>
#include <limits>
#include <vector>
#include "math_type.h"

namespace ornate {

template <typename T>
bool is_valid_price(T p) {
    return (p > 1e-6 && p < 1e8 && std::isfinite(p));
}

inline float calc_return(float p, float p0) {
    float r = NAN;
    if (is_valid_price(p) && is_valid_price(p0)) r = p / p0 - 1;
    if (!std::isfinite(r) || std::fabs(r) > 0.3) r = 0;
    return r * 10000;
}

inline bool is_valid_maturity(int n) { return (n > 0 && n < 9999); }
inline bool is_same_maturity(int mat0, int mat1) { return (mat0 == mat1 && is_valid_maturity(mat0)); }

template <typename T>
float calc_ir(const std::vector<T> &pnl, int32_t start_di = -1, int32_t end_di = -1) {
    if (pnl.empty()) return NAN;
    if (start_di >= end_di) return NAN;  // less than one day
    if (start_di < 0) start_di = 0;
    if (end_di < 0) end_di = static_cast<int32_t>(pnl.size() - 1);
    // delete leading NANs
    for (; start_di <= end_di; ++start_di) {
        if (pnl[start_di] != 0) break;
    }
    int32_t days = end_di - start_di + 1;
    if (days < 2) return NAN;
    float mean = 0;
    for (int32_t di = start_di; di <= end_di; ++di) mean += pnl[di];
    mean /= days;
    if (mean == 0) return 0;
    float var = 0;
    float t = 0;
    for (int32_t di = start_di; di <= end_di; ++di) {
        t = static_cast<float>(pnl[di] - mean);
        var += t * t;
    }
    var /= days;
    if (var <= 0) return NAN;
    return mean / sqrtf(var);
}

}  // namespace ornate

#endif
