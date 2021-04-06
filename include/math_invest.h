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

inline double get_valid_price_field(double val) {
    if (!ornate::is_valid_price(val))
        return NAN;
    else
        return val;
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
float calc_ir(const std::vector<T>& pnl, int32_t start_di = -1, int32_t end_di = -1) {
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

template <typename T>
inline double calc_accum_return(const std::vector<T>& rets) {
    double r = 1.0;
    for (const auto& ret : rets) {
        if (std::isfinite(ret)) {
            r *= (1.0 + ret);
        }
    }
    return r - 1.0;
}

/**
 * 十年收益 20%, total_ret = 0.2, n = 10
 */
inline double calc_avg_return(double total_ret, int n) {
    if (n <= 0)
        return 0;
    else
        return std::pow(1.0 + total_ret, 1.0 / n) - 1.0;
}

template <typename T>
inline std::vector<double> calc_bar_return_series(const std::vector<T>& signals, const std::vector<T>& rets,
                                                  int ins_num, bool is_signal_weighted, double open_t, double close_t) {
    double r = 1.0;
    int total = (int)signals.size();
    std::vector<int> ii2status(ins_num, 0);
    std::vector<double> ii2last_signal(ins_num, 0);
    std::vector<double> ret_vals;
    if (!is_signal_weighted) {
        std::fill(ii2last_signal.begin(), ii2last_signal.end(), 1.0 / ins_num);
    }
    for (int offset = 0; offset < total - ins_num; offset += ins_num) {
        double total_weight = 1;
        if (is_signal_weighted) {
            total_weight = 0;
            for (int ii = 0; ii < ins_num; ++ii) {
                if (isvalid(signals[offset + ii])) {
                    total_weight += std::abs(signals[offset + ii]);
                    ii2last_signal[ii] = signals[offset + ii];
                } else {
                    total_weight += std::abs(ii2last_signal[ii]);
                }
            }
        }

        if (total_weight > 1e-6) {
            for (int ii = 0; ii < ins_num; ++ii) {
                double ret = rets[offset + ins_num + ii];
                if (!std::isfinite(ret)) ret = 0;
                double sig = ii2last_signal[ii];

                if (sig > open_t || (ii2status[ii] > 0 && sig > close_t)) {
                    ii2status[ii] = 1;
                    r += (std::abs(sig) / total_weight) * (1. + ret);
                } else if (sig < -open_t || (ii2status[ii] < 0 && sig < -close_t)) {
                    ii2status[ii] = -1;
                    r += (std::abs(sig) / total_weight) * (1. - ret);
                } else {
                    ii2status[ii] = 0;
                }
            }
        }
        ret_vals.push_back(r);
    }
    return ret_vals;
}

template <typename T>
inline std::vector<std::vector<double>> calc_return_series_by_ii(const std::vector<T>& signals,
                                                                 const std::vector<T>& rets) {
    return {};
}

}  // namespace ornate

#endif
