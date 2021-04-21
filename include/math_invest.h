#ifndef ORNATE_MATH_INVEST_H
#define ORNATE_MATH_INVEST_H

#include <deque>
#include <limits>
#include <vector>
#include "math_type.h"
#include "math_vector.h"

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

template <typename T, typename T1>
std::vector<double> calc_bar_return_series(const T* signals, const T1* rets, int total, int ins_num,
                                           bool is_signal_weighted, double open_t, double close_t, int top_n = -1) {
    double nav = 1.0;
    std::vector<int> ii2status(ins_num, 0);
    std::vector<double> ret_vals;
    for (int offset = 0; offset < total; offset += ins_num) {
        std::vector<T> sigs_(signals + offset, signals + offset + ins_num);
        if (top_n > 0) {
            ornate::keep_top(sigs_, top_n, (T)0, true);
        }
        int valid_cnt = 0;
        double total_weight = 0;
        for (int ii = 0; ii < ins_num; ++ii) {
            if (isvalid(sigs_[ii])) {
                total_weight += std::abs(sigs_[ii]);
                ++valid_cnt;
            }
        }
        if (valid_cnt <= 0) continue;

        if (!is_signal_weighted) total_weight = 1.0;

        if (total_weight > 1e-6) {
            double tmp_nav = 0;
            for (int ii = 0; ii < ins_num; ++ii) {
                double ret = rets[offset + ii];
                if (!std::isfinite(ret)) ret = 0;
                double sig = sigs_[ii];
                if (!std::isfinite(sig)) sig = 0;
                double weight = std::abs(sig) / total_weight;
                if (!is_signal_weighted) weight = 1. / ins_num;

                if (sig > open_t || (ii2status[ii] > 0 && sig >= close_t)) {
                    ii2status[ii] = 1;
                    tmp_nav += nav * weight * (1. + ret);
                } else if (sig < -open_t || (ii2status[ii] < 0 && sig <= -close_t)) {
                    ii2status[ii] = -1;
                    tmp_nav += nav * weight * (1. - ret);
                } else {
                    ii2status[ii] = 0;
                    tmp_nav += nav * weight;
                }
            }
            nav = tmp_nav;
            if (nav < 0) nav = 0;
        }
        ret_vals.push_back(nav);
    }
    return ret_vals;
}

/**
 * @param signals
 * @param rets 未来一根bar收益
 * @param ins_num
 * @param is_signal_weighted
 * @param open_t open threshold
 * @param close_t close threshold
 * @return
 */
template <typename T, typename T1>
std::vector<double> calc_bar_return_series(const std::vector<T>& signals, const std::vector<T1>& rets, int ins_num,
                                           bool is_signal_weighted, double open_t, double close_t, int top_n = -1) {
    int total = (int)signals.size();
    return calc_bar_return_series(signals.data(), rets.data(), total, ins_num, is_signal_weighted, open_t, close_t,
                                  top_n);
}

template <typename T, typename T1>
std::vector<std::vector<double>> calc_return_series_by_ii(const std::vector<T>& signals, const std::vector<T1>& rets,
                                                          int ins_num, double open_t, double close_t) {
    std::vector<std::vector<double>> ret;
    for (int ii = 0; ii < ins_num; ++ii) {
        std::vector<T> ii_signal = ornate::skip_extract(signals, ins_num, ii);
        std::vector<T1> ii_ret = ornate::skip_extract(rets, ins_num, ii);
        ret.push_back(calc_bar_return_series(ii_signal, ii_ret, 1, false, open_t, close_t));
    }
    return ret;
}

template <typename T>
double calc_max_dropdown_ratio(const T* signals, int len) {
    double max_val = NAN, max_drop_ratio = NAN;
    for (int i = 0; i < len; ++i) {
        if (!std::isfinite(signals[i])) continue;
        if (std::isnan(max_val) || signals[i] > max_val) max_val = signals[i];
        double drop_ratio = 1. - signals[i] / max_val;
        if (std::isnan(max_drop_ratio) || drop_ratio > max_drop_ratio) max_drop_ratio = drop_ratio;
    }
    return max_drop_ratio;
}

template <typename T>
double calc_max_dropdown_ratio(const std::vector<T>& signals) {
    return calc_max_dropdown_ratio(signals.data(), signals.size());
}

}  // namespace ornate

#endif
