#pragma once

#include <cmath>

namespace ornate {
inline double Pow(double x, double y) {
  double s = (x > 1e-9) ? 1.0 : (x < -1e-9) ? -1.0 : 0.0;
  double ax = std::abs(x);
  if (y == 0.5) return s * std::sqrt(ax);
  if (y == 2.0) return s * ax * ax;
  if (y == 3.0) return s * ax * ax * ax;
  if (y == -1.0) return s / ax;
  return s * std::pow(ax, y);
}
inline double Sign(double x) {
    if (x > 1e-9) return 1.0;
    if (x < -1e-9) return -1.0;
    return 0.0;
}
inline double Na20(double x) { return !std::isfinite(x) ? 0.0 : x; }
inline double Inv(double x) { return 1.0 / x; }
inline double Sigmoid(double x) {
    if (x >= 36.0) return 1.0;
    if (x <= -36.0) return 0.0;
    return 1.0 / (1.0 + std::exp(-x));
}
inline double Relu(double x) { return x > 0 ? x : 0.0; }
inline double AbsSqrt(double x) {
    double s = Sign(x);
    return s * std::sqrt(std::abs(x));
}
}