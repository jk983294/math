#ifndef ORNATE_MATH_UTILS_H
#define ORNATE_MATH_UTILS_H

#include <deque>
#include <limits>
#include <string>
#include <vector>
#include "math_type.h"

namespace ornate {

constexpr double epsilon = 1e-9;

template <typename T>
bool IsValidData(T value) {
    return std::isfinite(value) && !std::isnan(value);
}

template <typename T1, typename T2>
bool FloatEqual(T1 a, T2 b) {
    if (IsValidData(a) && IsValidData(b)) {
        return std::fabs(a - b) < epsilon;
    } else if (std::isnan(a) && std::isnan(b)) {
        return true;
    } else if (!std::isfinite(a) && !std::isfinite(b)) {
        return true;
    } else {
        return false;
    }
}

template <typename T, typename T1>
void add_window_vector(std::vector<T>& y, size_t window, T1 val) {
    if (y.size() < window)
        y.push_back(val);
    else {
        for (size_t j = 1; j < window; ++j) {
            y[j - 1] = y[j];
        }
        y[window - 1] = val;
    }
}

/**
 * round_up(1.3456, 2) -> 1.35
 */
inline double round_up(double value, int decimal_places) {
    const double multiplier = std::pow(10.0, decimal_places);
    return static_cast<double>(std::lround(value * multiplier + .5)) / multiplier;
}

inline uint32_t NextPowerOf2(uint32_t n) {
    if (n && !(n & (n - 1))) return n;

    uint32_t count = 0;
    while (n != 0) {
        n >>= 1u;
        count += 1;
    }

    return 1u << count;
}

inline bool cmp_numeric(double l, double r) {
    if (std::isfinite(l) && std::isfinite(r))
        return l < r;
    else if (!std::isfinite(l) && !std::isfinite(r))
        return false;  // 永远让比较函数对相同元素返回false, otherwise violate Strict Weak Ordering
    else
        return !std::isfinite(l);
}

/**
 * hl越大，该alpha/decay系数越大，新值权重越小，老值权重越高，ma的时滞越大
 * hl越小，对新值越敏感
 * @param hl 半衰期
 */
inline double ema_hl2decay(double hl) { return pow(0.5, 1.0 / hl); }

inline bool is_same_sign(double val, int sign_) { return sign_ == 0 || sign_ * val > 0; }
inline bool is_same_sign(double val, double sign_) { return sign_ * val > 0; }
}  // namespace ornate

#endif
