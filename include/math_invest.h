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

/**
 * ii2status 1=多头, -1=空头, 0=空仓
 * sticky 可以通过设置 is_signal_weighted=false, open_t=close_t=0 来模拟
 * @param rets 未来一根bar的收益
 * @param total signals length
 * @param is_signal_weighted
 * @param open_t 多头开仓阈值, 空头开仓阈值=-open_t
 * @param close_t 多头平仓阈值, 空头平仓阈值=-close_t
 * @param top_n 每轮最大持有ins, -1表示没有限制
 * @param sticky true表示上一轮如果没有到平仓线,则保留,直到平仓阈值触发, false表示每轮只按signal强度,不考虑上轮持仓
 * @param cost 每笔的滑点佣金, 在每笔交易前收取
 */
template <typename T, typename T1>
std::vector<double> calc_bar_return_series(const T* signals, const T1* rets, int total, int ins_num,
                                           bool is_signal_weighted, double open_t, double close_t, int top_n = -1,
                                           bool sticky = true, double cost = 0) {
    double nav = 1.0;
    std::vector<int> ii2status(ins_num, 0);
    std::vector<double> ret_vals;
    for (int offset = 0; offset < total; offset += ins_num) {
        std::vector<T> sigs_(signals + offset, signals + offset + ins_num);
        int real_n = top_n;
        if (top_n > 0) {
            // real_n = ornate::keep_top(sigs_, top_n, (T)0, true);
            int tmp_n = 0;
            std::vector<std::pair<T, int>> sort_array;
            for (int ii = 0; ii < ins_num; ++ii) {
                if (isvalid(sigs_[ii])) {
                    // 没被选中, 但是没到平仓阈值,继续保留
                    if (sticky &&
                        ((ii2status[ii] > 0 && sigs_[ii] > close_t) || (ii2status[ii] < 0 && sigs_[ii] < -close_t))) {
                        ++tmp_n;
                        continue;
                    }
                    sort_array.emplace_back(std::abs(sigs_[ii]), ii);
                }
            }

            int needed = top_n - tmp_n;
            if (needed < 0) {
                printf("should not happen");
                needed = 0;
            }
            if (needed < (int)sort_array.size()) {
                std::sort(sort_array.begin(), sort_array.end(),
                          [](const auto& l, const auto& r) { return l.first > r.first; });
                for (int i = needed; i < (int)sort_array.size(); ++i) {
                    sigs_[sort_array[i].second] = 0;
                }
                real_n = top_n;
            } else {
                real_n = tmp_n + (int)sort_array.size();
            }
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
            int selected_cnt = 0;
            for (int ii = 0; ii < ins_num; ++ii) {
                double ret = rets[offset + ii];
                if (!std::isfinite(ret)) ret = 0;
                double sig = sigs_[ii];
                if (!std::isfinite(sig)) sig = 0;
                double weight = std::abs(sig) / total_weight;
                if (!is_signal_weighted) {
                    if (real_n > 0)
                        weight = 1. / real_n;
                    else
                        weight = 1. / ins_num;
                }

                if (sig > open_t || (sticky && ii2status[ii] > 0 && sig > close_t)) {
                    double tmp_cost = 0;
                    if (ii2status[ii] <= 0) tmp_cost = cost;  // 老仓位方向不同
                    ii2status[ii] = 1;
                    tmp_nav += nav * weight * (1. + ret - tmp_cost);
                    ++selected_cnt;
                } else if (sig < -open_t || (sticky && ii2status[ii] < 0 && sig < -close_t)) {
                    double tmp_cost = 0;
                    if (ii2status[ii] >= 0) tmp_cost = cost;  // 老仓位方向不同
                    ii2status[ii] = -1;
                    tmp_nav += nav * weight * (1. - ret - tmp_cost);
                    ++selected_cnt;
                } else {
                    ii2status[ii] = 0;
                    if (!is_signal_weighted && real_n > 0 && sig == 0)
                        tmp_nav += 0;
                    else {
                        tmp_nav += nav * weight;
                        ++selected_cnt;
                    }
                }
            }
            if (!is_signal_weighted && real_n > 0) {
                if (selected_cnt < real_n) {
                    tmp_nav += nav * (real_n - selected_cnt) / real_n;
                } else if (selected_cnt > real_n) {
                    printf("error, should not happen");
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
                                           bool is_signal_weighted, double open_t, double close_t, int top_n = -1,
                                           bool sticky = true, double cost = 0) {
    int total = (int)signals.size();
    return calc_bar_return_series(signals.data(), rets.data(), total, ins_num, is_signal_weighted, open_t, close_t,
                                  top_n, sticky, cost);
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
