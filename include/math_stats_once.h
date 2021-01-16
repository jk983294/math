#ifndef ORNATE_MATH_STATS_ONCE_H
#define ORNATE_MATH_STATS_ONCE_H

#include "math_utils.h"

namespace ornate {
struct rolling_mean_once {
    double total_sum{0};
    int cnt{0};

    void operator()(double x) {
        if (std::isfinite(x)) {
            total_sum += x;
            ++cnt;
        }
    }

    void clear() {
        total_sum = 0;
        cnt = 0;
    }

    double final() {
        if (cnt > 0)
            return total_sum / cnt;
        else
            return NAN;
    }
};

struct rolling_sd_once {
    double total_sum{0}, total_square_sum{0};
    int cnt{0};

    void operator()(double x) {
        if (std::isfinite(x)) {
            total_sum += x;
            total_square_sum += x * x;
            ++cnt;
        }
    }

    void clear() {
        total_sum = 0;
        cnt = 0;
    }

    double final() {
        if (cnt > 1) {
            double mean = total_sum / cnt;
            double variance = (total_square_sum - mean * mean * cnt) / (cnt - 1);
            return std::sqrt(variance);
        } else {
            return NAN;
        }
    }
};
}  // namespace ornate

#endif
