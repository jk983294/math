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
 * @param pre_close_ret_t 提前close仓位的条件,如果收益大于该阈值
 * @param pre_close_tick_t 提前close仓位的条件,如果持仓时间大于该阈值
 */
struct TickSimStat {
    TickSimStat(int total, int ins_num) : m_total{total}, m_ins_num{ins_num} {
        ii2status.resize(m_ins_num, 0);
        hold_ret.resize(m_ins_num, 0);
        sigs_.resize(m_ins_num, NAN);
        ret_vals.reserve(m_total / m_ins_num);
        hold_ticks.reserve(m_total);
        m_rets.reserve(m_total);
    }

    template <typename T, typename T1>
    std::vector<double> calc_bar_return_series(const T* signals, const T1* rets) {
        double nav = 1.0;
        std::fill(ii2status.begin(), ii2status.end(), 0);
        ret_vals.clear();
        hold_ticks.clear();
        m_rets.clear();
        for (int offset = 0, ti = 0; offset < m_total; offset += m_ins_num, ++ti) {
            if (m_ti_num > 0 && (ti + 1) % m_ti_num == 0) ti = 0;
            bool clear_ = m_ti_num > 0 && !m_clear_ticks.empty() && m_clear_ticks[ti];
            if (clear_) {
                std::fill(sigs_.begin(), sigs_.end(), NAN);
            } else {
                std::copy(signals + offset, signals + offset + m_ins_num, sigs_.begin());
                if (m_is_neutral) {
                    double mean_ = mean(sigs_);
                    vs_minus_inplace(sigs_, mean_);
                }
            }
            int real_n = m_top_n;
            if (m_top_n > 0) {
                real_n = choose_top_n();
            }
            double total_weight = calc_weight();

            if (total_weight > 1e-6) {
                double tmp_nav = 0;
                int selected_cnt = 0;
                for (int ii = 0; ii < m_ins_num; ++ii) {
                    double ret = rets[offset + ii];
                    if (!std::isfinite(ret)) ret = 0;
                    double sig = sigs_[ii];
                    if (!std::isfinite(sig)) sig = 0;
                    double weight = std::abs(sig) / total_weight;
                    if (!m_is_signal_weighted) {
                        if (real_n > 0)
                            weight = 1. / real_n;
                        else {
                            weight = 1. / m_ins_num;  // top_n = -1, 所有ii等权重
                        }
                    }

                    if (sig > m_open_t || (m_sticky && ii2status[ii] > 0 && sig > m_close_t)) {
                        double tmp_cost = 0;
                        if (ii2status[ii] <= 0) tmp_cost = m_cost;  // 老仓位方向不同
                        if (ii2status[ii] < 0) {                    // 老仓位方向不同
                            record_ret_and_ticks(ii);
                        }
                        ii2status[ii] += 1;
                        accum_ret(ii, ret);
                        tmp_nav += nav * weight * (1. + ret - tmp_cost);
                        ++selected_cnt;
                    } else if (sig < -m_open_t || (m_sticky && ii2status[ii] < 0 && sig < -m_close_t)) {
                        double tmp_cost = 0;
                        if (ii2status[ii] >= 0) tmp_cost = m_cost;  // 老仓位方向不同
                        if (ii2status[ii] > 0) {                    // 老仓位方向不同
                            record_ret_and_ticks(ii);
                        }
                        ii2status[ii] -= 1;
                        accum_ret(ii, ret);
                        tmp_nav += nav * weight * (1. - ret - tmp_cost);
                        ++selected_cnt;
                    } else {
                        record_ret_and_ticks(ii);
                        if (!m_is_signal_weighted && real_n > 0 && sig == 0)
                            tmp_nav += 0;
                        else {
                            tmp_nav += nav * weight;
                            ++selected_cnt;
                        }
                    }
                }
                if (!m_is_signal_weighted && real_n > 0) {
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
            if (clear_) clear_pos();
        }

        clear_pos();
        return ret_vals;
    }

    std::vector<double> calc_navs(const std::vector<double>& signals, const std::vector<double>& rets) {
        return calc_bar_return_series(signals.data(), rets.data());
    }

    void set_is_signal_weighted(bool is_signal_weighted) { m_is_signal_weighted = is_signal_weighted; }
    void set_open_t(double open_t) { m_open_t = open_t; }
    void set_close_t(double close_t) { m_close_t = close_t; }
    void set_cost(double cost) { m_cost = cost; }
    void set_pre_close_ret_t(double pre_close_ret_t) { m_pre_close_ret_t = pre_close_ret_t; }
    void set_top_n(int top_n) { m_top_n = top_n; }
    void set_pre_close_tick_t(int pre_close_tick_t) { m_pre_close_tick_t = pre_close_tick_t; }
    void set_sticky(bool sticky) { m_sticky = sticky; }
    void set_is_neutral(bool is_neutral) {
        m_is_neutral = is_neutral;
        if (m_is_neutral && m_top_n % 2 == 1) m_top_n -= 1;
    }
    void set_ti_num(int ti_num) {
        m_ti_num = ti_num;
        if (m_ti_num > 0) m_clear_ticks.resize(m_ti_num, false);
    }
    const std::vector<int>& get_hold_ticks() { return hold_ticks; }
    const std::vector<double>& get_rets() { return m_rets; }
    const std::vector<double>& get_navs() { return ret_vals; }
    void set_clear_ticks(const std::vector<int>& ticks) {
        if (m_ti_num > 0) {
            m_clear_ticks.resize(m_ti_num);
            std::fill(m_clear_ticks.begin(), m_clear_ticks.end(), false);
            for (int ti : ticks) {
                if (ti < 0 || ti >= m_ti_num) continue;
                m_clear_ticks[ti] = true;
            }
        }
    }

private:
    void accum_ret(int ii, double ret_) {
        if (std::isfinite(ret_)) {
            if (ii2status[ii] < 0) {
                m_rets[ii] -= ret_;
            } else if (ii2status[ii] > 0) {
                m_rets[ii] += ret_;
            }
        }
    }
    void record_ret_and_ticks(int ii) {
        if (ii2status[ii] == 0) {
            return;
        } else if (ii2status[ii] < 0) {
            hold_ticks.push_back(-ii2status[ii]);
            m_rets.push_back(hold_ret[ii]);
        } else if (ii2status[ii] > 0) {
            hold_ticks.push_back(ii2status[ii]);
            m_rets.push_back(hold_ret[ii]);
        }
        ii2status[ii] = 0;
        hold_ret[ii] = 0;
    }
    int choose_top_n() {
        if (m_is_neutral) return neutral_choose_top_n();
        int tmp_n = 0;
        std::vector<std::pair<double, int>> sort_array;
        for (int ii = 0; ii < m_ins_num; ++ii) {
            if (isvalid(sigs_[ii])) {
                // 没被选中, 但是没到平仓阈值,继续保留
                if (m_sticky &&
                    ((ii2status[ii] > 0 && sigs_[ii] > m_close_t) || (ii2status[ii] < 0 && sigs_[ii] < -m_close_t))) {
                    if (not((std::isfinite(m_pre_close_ret_t) && m_rets[ii] > m_pre_close_ret_t) ||
                            (m_pre_close_tick_t > 0 && std::abs(ii2status[ii]) > m_pre_close_tick_t))) {
                        ++tmp_n;
                        continue;
                    }
                }
                sort_array.emplace_back(std::abs(sigs_[ii]), ii);
            }
        }

        int needed = m_top_n - tmp_n;  // 除了需要保留的，还需要多少补充
        if (needed < 0) {
            printf("should not happen");
            needed = 0;
        }
        if (needed < (int)sort_array.size()) {
            std::sort(sort_array.begin(), sort_array.end(),
                      [](const auto& l, const auto& r) { return l.first > r.first; });
            for (int i = needed; i < (int)sort_array.size(); ++i) {
                sigs_[sort_array[i].second] = NAN;  // 消除掉不需要 ii signal
            }
            return m_top_n;
        } else {
            return tmp_n + (int)sort_array.size();
        }
    }
    int neutral_choose_top_n() {
        int pos_n = 0, neg_n = 0;
        std::vector<std::pair<double, int>> sort_array_pos;
        std::vector<std::pair<double, int>> sort_array_neg;
        for (int ii = 0; ii < m_ins_num; ++ii) {
            if (isvalid(sigs_[ii])) {
                // 没被选中, 但是没到平仓阈值,继续保留
                if (m_sticky) {
                    if (ii2status[ii] > 0 && sigs_[ii] > m_close_t) {
                        ++pos_n;
                    } else if (ii2status[ii] < 0 && sigs_[ii] < -m_close_t) {
                        ++neg_n;
                    } else if (sigs_[ii] > 0.) {
                        sort_array_pos.emplace_back(std::abs(sigs_[ii]), ii);
                    } else {
                        sort_array_neg.emplace_back(std::abs(sigs_[ii]), ii);
                    }
                }
            }
        }

        std::sort(sort_array_pos.begin(), sort_array_pos.end(),
                  [](const auto& l, const auto& r) { return l.first > r.first; });
        int pos_need = m_top_n / 2 - pos_n;
        if (pos_need < (int)sort_array_pos.size()) {
            for (int i = pos_need; i < (int)sort_array_pos.size(); ++i) {
                sigs_[sort_array_pos[i].second] = NAN;  // 消除掉不需要 ii signal
            }
            pos_n += pos_need;
        } else {
            pos_n += (int)sort_array_pos.size();
        }

        std::sort(sort_array_neg.begin(), sort_array_neg.end(),
                  [](const auto& l, const auto& r) { return l.first > r.first; });
        int neg_need = m_top_n / 2 - pos_n;
        if (neg_need < (int)sort_array_neg.size()) {
            for (int i = neg_need; i < (int)sort_array_neg.size(); ++i) {
                sigs_[sort_array_neg[i].second] = NAN;  // 消除掉不需要 ii signal
            }
            pos_n += neg_need;
        } else {
            pos_n += (int)sort_array_neg.size();
        }
        return pos_n + neg_n;
    }

    double calc_weight() {
        int valid_cnt = 0;
        double total_weight = 0;
        for (int ii = 0; ii < m_ins_num; ++ii) {
            if (isvalid(sigs_[ii])) {
                total_weight += std::abs(sigs_[ii]);
                ++valid_cnt;
            }
        }
        if (valid_cnt <= 0) {
            return 0;
        } else {
            if (!m_is_signal_weighted)
                return 1.0;
            else
                return total_weight;
        }
    }

    void clear_pos() {
        for (int ii = 0; ii < m_ins_num; ++ii) {
            record_ret_and_ticks(ii);
        }
    }

private:
    int m_total{0};
    int m_ins_num{0};
    int m_ti_num{-1};
    int m_top_n{-1};
    int m_pre_close_tick_t{-1};
    bool m_is_signal_weighted{false};
    bool m_sticky{true};
    bool m_is_neutral{false};
    double m_open_t{0};
    double m_close_t{0};
    double m_cost{0};
    double m_pre_close_ret_t{NAN};
    std::vector<int> ii2status, hold_ticks;
    std::vector<bool> m_clear_ticks;
    std::vector<double> ret_vals;
    std::vector<double> sigs_;
    std::vector<double> hold_ret, m_rets;
};

template <typename T, typename T1>
std::vector<double> calc_bar_return_series(const T* signals, const T1* rets, int total, int ins_num,
                                           bool is_signal_weighted, double open_t, double close_t, int top_n = -1,
                                           bool sticky = true, double cost = 0, bool neutral = false) {
    TickSimStat stat(total, ins_num);
    stat.set_is_signal_weighted(is_signal_weighted);
    stat.set_open_t(open_t);
    stat.set_close_t(close_t);
    stat.set_top_n(top_n);
    stat.set_sticky(sticky);
    stat.set_cost(cost);
    stat.set_is_neutral(neutral);
    stat.calc_bar_return_series(signals, rets);
    return stat.get_navs();
}

template <typename T, typename T1>
std::vector<int> calc_hold_tick(const T* signals, const T1* rets, int total, int ins_num, double open_t, double close_t,
                                int top_n = -1, bool sticky = true, bool neutral = false) {
    TickSimStat stat(total, ins_num);
    stat.set_open_t(open_t);
    stat.set_close_t(close_t);
    stat.set_top_n(top_n);
    stat.set_sticky(sticky);
    stat.set_is_neutral(neutral);
    stat.calc_bar_return_series(signals, rets);
    return stat.get_hold_ticks();
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
                                           bool sticky = true, double cost = 0, bool neutral = false) {
    int total = (int)signals.size();
    return calc_bar_return_series(signals.data(), rets.data(), total, ins_num, is_signal_weighted, open_t, close_t,
                                  top_n, sticky, cost, neutral);
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
std::pair<double, double> calc_max_dropdown_ratio(const T* signals, int len) {
    double max_val = NAN, max_drop_ratio = NAN, time_ratio = NAN;
    int max_idx = -1;
    for (int i = 0; i < len; ++i) {
        if (!std::isfinite(signals[i])) continue;
        if (std::isnan(max_val) || signals[i] > max_val) {
            max_val = signals[i];
            max_idx = i;
        }
        double drop_ratio = 1. - signals[i] / max_val;
        if (std::isnan(max_drop_ratio) || drop_ratio > max_drop_ratio) {
            max_drop_ratio = drop_ratio;
            time_ratio = double(i - max_idx + 1) / len;
        }
    }
    return {max_drop_ratio, time_ratio};
}

template <typename T>
std::pair<double, double> calc_max_dropdown_ratio(const std::vector<T>& signals) {
    return calc_max_dropdown_ratio(signals.data(), signals.size());
}

/**
 * m_status 1=多头, -1=空头, 0=空仓
 * @param rets 未来一根bar的收益
 * @param open_t 多头开仓阈值, 空头开仓阈值=-open_t
 * @param close_t 多头平仓阈值, 空头平仓阈值=-close_t
 * @param cost 每笔的滑点佣金, 在每笔交易前收取
 * @param m_profit_ratio 提前close仓位的条件,如果收益大于该阈值
 * @param m_stop_ratio 提前close仓位的条件,如果 loss > 该阈值
 * @param m_stop_tick 提前close仓位的条件,如果持仓时间大于该阈值
 */
struct TickSimStatSingle {
    TickSimStatSingle() = default;
    TickSimStatSingle(double open_t, double close_t, double cost, double stop_ratio, double profit_ratio, int stop_tick)
        : m_open_t{open_t},
          m_close_t{close_t},
          m_cost{cost},
          m_stop_ratio{stop_ratio},
          m_profit_ratio{profit_ratio},
          m_stop_tick{stop_tick} {}

    std::vector<double> calc_navs(const std::vector<double>& signals, const std::vector<double>& rets,
                                  double init_nav = 1.0) {
        return calc_bar_return_series(signals.data(), rets.data(), signals.size(), init_nav);
    }

    template <typename T, typename T1>
    std::vector<double> calc_bar_return_series(const T* signals, const T1* rets, int len, double init_nav = 1.0) {
        clear(len);
        m_nav = init_nav;
        double curr_nav = m_nav;
        accum_ret = max_ret = min_ret = 0.0;
        m_status = 0;
        m_pre_stop_dir = 0;
        for (int i = 0; i < len; ++i) {
            double next_ret1 = (std::isfinite(rets[i]) ? rets[i] : 0.);
            int curr_dir = get_signal(signals[i]);

            if (has_pos(m_status)) {
                if (has_long_pos(m_status)) {
                    curr_nav = m_nav * (1.0 + accum_ret - m_cost);
                } else if (has_short_pos(m_status)) {
                    curr_nav = m_nav * (1.0 - accum_ret - m_cost);
                }

                if (is_opposite_sig(m_status, curr_dir)) {
                    clear_pos();
                    open_pos(curr_dir, next_ret1);
                } else {  // check drop ratio
                    if (is_stop_win() || is_stop_loss() || is_stop_time()) {
                        m_pre_stop_dir = m_status;
                        clear_pos();
                    } else if (m_status == 0) {
                        clear_pos();
                    } else {
                        m_status += curr_dir;
                        accum_ret = (1.0 + accum_ret) * (1.0 + next_ret1) - 1.0;
                        max_ret = std::max(max_ret, accum_ret);
                        min_ret = std::min(min_ret, accum_ret);
                    }
                }
            } else if ((is_long_signal(curr_dir) && !is_long_signal(m_pre_stop_dir)) ||
                       (is_short_signal(curr_dir) && !is_short_signal(m_pre_stop_dir))) {
                open_pos(curr_dir, next_ret1);
            }
            m_navs[i] = curr_nav;
        }

        clear_pos();
        return m_navs;
    }

    void set_open_t(double open_t) { m_open_t = open_t; }
    void set_close_t(double close_t) { m_close_t = close_t; }
    void set_cost(double cost) { m_cost = cost; }
    void set_stop_ratio(double stop_ratio) { m_stop_ratio = stop_ratio; }
    void set_profit_ratio(double profit_ratio) { m_profit_ratio = profit_ratio; }
    void set_stop_tick(int stop_tick) { m_stop_tick = stop_tick; }
    const std::vector<int>& get_hold_ticks() { return m_hold_ticks; }
    const std::vector<double>& get_navs() { return m_navs; }

private:
    void clear(int len) {
        m_navs.clear();
        m_navs.resize(len, 0);
        m_hold_ticks.clear();
        m_hold_ticks.reserve(len);
    }

    static bool has_long_pos(int pos) { return pos > 0; }
    static bool has_short_pos(int pos) { return pos < 0; }
    static bool has_pos(int pos) { return has_long_pos(pos) || has_short_pos(pos); }
    static bool is_long_signal(int pos) { return pos > 0; }
    static bool is_short_signal(int pos) { return pos < 0; }
    int get_signal(double sig) const {
        if (!std::isfinite(sig)) {  // if NAN signal, return last signal
            if (m_status > 0)
                return 1;
            else if (sig < 0)
                return -1;
            return 0;
        } else {
            if (sig > m_open_t || (m_status > 0 && sig > m_close_t))
                return 1;
            else if (sig < -m_open_t || (m_status < 0 && sig < -m_close_t))
                return -1;
            return 0;
        }
    }
    static bool is_opposite_sig(int old_sig, int sig) { return old_sig * sig < 0; }
    bool is_stop_win() const {
        if (std::isfinite(m_profit_ratio) && m_profit_ratio > 1e-6) {
            if (has_long_pos(m_status) && accum_ret > m_profit_ratio) {
                return true;
            } else if (has_short_pos(m_status) && -accum_ret < -m_profit_ratio) {
                return true;
            }
        }
        return false;
    }
    bool is_stop_loss() const {
        if (std::isfinite(m_stop_ratio) && m_stop_ratio > 1e-6) {
            if (has_long_pos(m_status)) {
                return max_ret - accum_ret > m_stop_ratio;
            } else if (has_short_pos(m_status)) {
                return accum_ret - min_ret > m_stop_ratio;
            }
        }
        return false;
    }

    bool is_stop_time() const { return m_stop_tick > 0 && std::abs(m_status) >= m_stop_tick; }

    void clear_pos() {
        if (has_long_pos(m_status)) {
            m_nav *= (1.0 + accum_ret - m_cost);
            m_hold_ticks.push_back(m_status);
        } else if (has_short_pos(m_status)) {
            m_nav *= (1.0 - accum_ret - m_cost);
            m_hold_ticks.push_back(-m_status);
        }
        m_status = 0;
        accum_ret = max_ret = min_ret = 0.0;
    }

    void open_pos(int curr_pos, double curr_ret) {
        m_status = curr_pos;
        m_pre_stop_dir = 0;
        accum_ret = curr_ret;
        max_ret = std::max(curr_ret, 0.);
        min_ret = std::min(curr_ret, 0.);
    }

private:
    double m_open_t{0};
    double m_close_t{0};
    double m_cost{0};
    double m_stop_ratio{NAN};
    double m_profit_ratio{NAN};
    double accum_ret = 0, max_ret = 0, min_ret = 0;
    double m_nav{0};
    int m_stop_tick{-1};
    int m_status{0};
    int m_pre_stop_dir{0};  // if stop, then succeeding same dir signal will be ignored
    std::vector<int> m_hold_ticks;
    std::vector<double> m_navs;
};

template <typename T, typename T1>
std::vector<double> hft_calc_bar_return_series(const T* signals, const T1* rets, int len, double open_t, double close_t,
                                               double cost = 0, double stop_ratio = NAN, double profit_ratio = NAN,
                                               int stop_tick = -1) {
    TickSimStatSingle stat(open_t, close_t, cost, stop_ratio, profit_ratio, stop_tick);
    return stat.calc_bar_return_series(signals, rets, len);
}

template <typename T, typename T1>
std::vector<double> hft_calc_bar_return_series_vec(const std::vector<T>& signals, const std::vector<T1>& rets,
                                                   double open_t, double close_t, double cost = 0,
                                                   double stop_ratio = NAN, double profit_ratio = NAN,
                                                   int stop_tick = -1) {
    return hft_calc_bar_return_series(signals.data(), rets.data(), signals.size(), open_t, close_t, cost, stop_ratio,
                                      profit_ratio, stop_tick);
}

inline bool sign_ratio_out_range(double neg_ratio, double sign_threshold) {
    return std::isfinite(sign_threshold) &&
           (std::isnan(neg_ratio) || neg_ratio < std::min(sign_threshold, 1.0 - sign_threshold) ||
            neg_ratio > std::max(sign_threshold, 1.0 - sign_threshold));
}

}  // namespace ornate

#endif
