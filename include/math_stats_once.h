#ifndef ORNATE_MATH_STATS_ONCE_H
#define ORNATE_MATH_STATS_ONCE_H

#include <math_stats.h>
#include <math_utils.h>

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

struct rolling_sm_once {
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

    std::pair<double, double> final() {
        double mean = NAN, sd = NAN;
        if (cnt > 0) {
            mean = total_sum / cnt;
        }
        if (cnt > 1) {
            double variance = (total_square_sum - mean * mean * cnt) / (cnt - 1);
            sd = std::sqrt(variance);
        }
        return {mean, sd};
    }
};

struct rolling_hl_once {
    double high{NAN}, low{NAN};

    void operator()(double x) {
        if (std::isfinite(x)) {
            if (std::isnan(high)) {
                high = x;
                low = x;
            } else if (x < low) {
                low = x;
            } else if (x > high) {
                high = x;
            }
        }
    }

    void clear() {
        high = NAN;
        low = NAN;
    }

    std::pair<double, double> final() { return {high, low}; }
};

struct regression2_once {
    long double sum_x{0}, sum_y{0}, sum_x2{0}, sum_xy{0};
    int m_valid_count{0};
    double a{NAN}, b{NAN}, mean_y, res_squared{0}, y_diff_squared{0};

    explicit regression2_once() {}

    void clear() {
        a = b = NAN;
        sum_x = sum_y = sum_x2 = sum_xy = 0;
        m_valid_count = 0;
    }

    void calc_coef() {
        a = b = NAN;
        if (m_valid_count > 1) {
            b = (m_valid_count * sum_xy - sum_x * sum_y) / (m_valid_count * sum_x2 - sum_x * sum_x);
            a = (sum_y - b * sum_x) / m_valid_count;
        }
    }

    double calc_fitted(double x) const { return a + b * x; }
    double calc_residual(double x, double y) const { return y - calc_fitted(x); }
    void r2_pre_calc() {
        calc_coef();
        res_squared = y_diff_squared = 0;
        if (m_valid_count > 1) {
            mean_y = sum_y / m_valid_count;
        } else {
            mean_y = NAN;
        }
    }
    void r2_single(double x, double y) {
        if (m_valid_count <= 1) return;
        if (!std::isfinite(x) || !std::isfinite(y)) return;
        double res = calc_residual(x, y);
        res_squared += res * res;
        double diff = y - mean_y;
        y_diff_squared += diff * diff;
    }
    double get_r2() const {
        if (m_valid_count > 1) {
            return 1. - res_squared / y_diff_squared;
        } else {
            return NAN;
        }
    }

    void operator()(const double y, const double x) {
        if (std::isfinite(y) && std::isfinite(x)) {
            sum_x += x;
            sum_y += y;
            sum_x2 += x * x;
            sum_xy += x * y;
            ++m_valid_count;
        }
    }

    void init() { clear(); }
};

template <typename T>
std::vector<double> lr_fitted(const T* y_vec, const T* x_vec, int num) {
    std::vector<double> pred(num, NAN);
    ornate::regression2_once f;
    double mean_ = ornate::mean(x_vec, num);
    double sd_ = ornate::std(x_vec, num);
    for (int i = 0; i < num; ++i) {
        if (!ornate::isvalid(x_vec[i]) || !ornate::isvalid(y_vec[i])) {
            continue;
        }
        f(y_vec[i], cap_within_sd(x_vec[i], mean_, sd_));
    }

    f.calc_coef();

    for (int i = 0; i < num; ++i) {
        if (ornate::isvalid(x_vec[i])) {
            pred[i] = f.calc_fitted(x_vec[i]);
        }
    }
    return pred;
}

}  // namespace ornate

#endif
