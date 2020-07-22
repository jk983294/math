#ifndef ORNATE_MATH_UTILS_H
#define ORNATE_MATH_UTILS_H

#include <deque>
#include <limits>
#include <string>
#include <vector>
#include "math_type.h"

namespace ornate {

constexpr double epsilon = 1e-6;

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
    } else {
        return false;
    }
}

inline double round_up(double value, int decimal_places) {
    const double multiplier = std::pow(10.0, decimal_places);
    return static_cast<double>(std::lround(value * multiplier + .5)) / multiplier;
}
}  // namespace ornate

#endif
