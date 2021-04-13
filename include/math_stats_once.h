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
        total_square_sum = total_sum = 0;
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

struct rolling_pcor_once {
    double sumx{0}, sumy{0}, sumxy{0}, sum_x2{0}, sum_y2{0};
    int m_valid_count{0};

    rolling_pcor_once() = default;

    void clear() {
        sumx = sumy = sumxy = sum_x2 = sum_y2 = 0;
        m_valid_count = 0;
    }

    void add_new(double data0, double data1) {
        if (std::isfinite(data0) && std::isfinite(data1)) {
            sumxy += data0 * data1;
            sumx += data0;
            sumy += data1;
            sum_x2 += data0 * data0;
            sum_y2 += data1 * data1;
            ++m_valid_count;
        }
    }

    double final() {
        if (m_valid_count > 1) {
            double mean_x = sumx / m_valid_count;
            double mean_y = sumy / m_valid_count;
            double cov = (sumxy - mean_x * mean_y * m_valid_count) / (m_valid_count - 1);
            double var1 = (sum_x2 - mean_x * mean_x * m_valid_count) / (m_valid_count - 1);
            double var2 = (sum_y2 - mean_y * mean_y * m_valid_count) / (m_valid_count - 1);
            if (std::isfinite(var1) && std::isfinite(var2)) {
                double numerator = sqrt(var1) * sqrt(var2);
                if (numerator >= epsilon) return cov / numerator;
            }
        }
        return NAN;
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

struct rolling_all_once {
    long double total_sum{0}, total_square_sum{0}, total_cube_sum{0};
    double high{NAN}, low{NAN}, last{NAN}, first{NAN};
    int cnt{0};

    void operator()(double x) {
        if (std::isfinite(x)) {
            last = x;
            if (std::isnan(high)) {
                high = x;
                low = x;
            } else if (x < low) {
                low = x;
            } else if (x > high) {
                high = x;
            }

            total_sum += x;
            total_square_sum += x * x;
            total_cube_sum += x * x * x;
            ++cnt;

            if (cnt == 1) {
                first = x;
            }
        }
    }

    void clear() {
        high = NAN;
        low = NAN;
        last = NAN;
        total_cube_sum = total_square_sum = total_sum = 0;
        cnt = 0;
    }

    std::pair<double, double> get_high_low() { return {high, low}; }

    std::pair<double, double> get_mean_sd() {
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

    double get_rsd() {
        double mean, sd;
        std::tie(mean, sd) = get_mean_sd();
        if (std::abs(mean) > 1e-9)
            return sd / mean;
        else
            return NAN;
    }

    double get_skew() {
        if (cnt >= 2) {
            double mean = total_sum / cnt;
            double mean2 = mean * mean;
            double var = total_square_sum / cnt - mean2;
            if (var <= 1e-14)
                return NAN;
            else {
                double mean3 = mean * mean * mean;
                double m3 = total_cube_sum / cnt - 3 * mean * total_square_sum / cnt + 2 * mean3;
                return m3 / std::pow(var, 1.5);
            }
        } else
            return NAN;
    }
};

template <typename T>
std::vector<double> lr_fitted(const T* y_vec, const T* x_vec, int num, bool is_cap = true) {
    std::vector<double> pred(num, NAN);
    ornate::regression2_once f;
    double mean_ = ornate::mean(x_vec, num);
    double sd_ = ornate::std(x_vec, num);
    for (int i = 0; i < num; ++i) {
        if (!ornate::isvalid(x_vec[i]) || !ornate::isvalid(y_vec[i])) {
            continue;
        }
        double x_ = x_vec[i];
        if (is_cap) x_ = cap_within_sd(x_vec[i], mean_, sd_);
        f(y_vec[i], x_);
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
