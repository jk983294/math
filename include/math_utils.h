#ifndef ORNATE_MATH_UTILS_H
#define ORNATE_MATH_UTILS_H

#include <deque>
#include <limits>
#include <vector>
#include "math_type.h"

namespace ornate {

constexpr double epsilon = 1e-6;

template <typename T>
bool IsValidData(T value) {
    return std::isfinite(value) && !std::isnan(value);
}

template <typename T>
bool FloatEqual(T a, T b) {
    if (IsValidData(a) && IsValidData(b)) {
        return std::fabs(a - b) < epsilon;
    } else if (std::isnan(a) && std::isnan(b)) {
        return true;
    } else {
        return false;
    }
}
}  // namespace ornate

#endif
