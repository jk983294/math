#include <math_dummy.h>
#include <math_vector.h>
#include <algorithm>
#include <cmath>
#include <cstddef>

namespace ornate {
int dummy_add(int x, int y) { return x + y; }
int universal_answer() { return 42; }
double dummy_ts_cross(const double* x_, const double* y, std::size_t i, std::size_t n) {
    const std::size_t size = i + 1;
    if (size < n + 2) {
        return NAN;
    }
    double x0 = x_[i - n - 1];
    double y0 = y[i - n - 1];
    if (std::isnan(x0) || std::isnan(y0)) {
        return NAN;
    } else {
        double x1_ = x_[i - n];
        double y1 = y[i - n];
        if (std::isnan(x1_) || std::isnan(y1)) {
            return NAN;
        } else if (x0 < y0 && x1_ > y1) {
            return 1;
        } else if (x0 > y0 && x1_ < y1) {
            return -1;
        } else {
            return 0;
        }
    }
}

double dummy_ts_argmax(const double* x, std::size_t i, std::size_t n) {
    const std::size_t size = i + 1;
    const std::size_t N = (n == 0) ? size : std::min(size, n);
    int max_idx = -1;
    double max_val = NAN;
    for (std::size_t j = 0; j < N; j++) {
        double curr_val = x[i - N + j + 1];
        if (std::isfinite(curr_val) && (std::isnan(max_val) || max_val < curr_val)) {
            max_val = curr_val;
            max_idx = int(j);
        }
    }
    if (N == 1) {
        if (max_idx >= 0) {
            return 1.;
        } else
            return NAN;
    } else {
        if (max_idx >= 0) {
            return double(N - 1 - max_idx) / double(N - 1);
        } else
            return NAN;
    }
}

double dummy_ts_backward_cpn(const double* x, std::size_t i, std::size_t n, int sign) {
    const std::size_t size = i + 1;
    const std::size_t N = (n == 0) ? size : std::min(size, n);
    double num = 0;
    for (std::size_t j = 0; j < N; j++) {
        if (std::isnan(x[i - j]) || (sign == 0 ? x[i] : sign) * x[i - j] <= 0) {
            break;
        } else {
            num++;
        }
    }
    return num;
}

double dummy_ts_acp(const double* x, int i, int n, int lag) {
    int bg = i - n + 1;
    int end = i + 1;

    bg = std::max(bg, 0);
    std::size_t valid_cnt = 0;
    std::size_t positive_cnt = 0;
    for (int j = bg; j < end - lag; j++) {
        const double vx = x[j];
        const double vx2 = x[j + lag];
        if (std::isnan(vx) or std::isnan(vx2)) continue;
        valid_cnt++;
        double res = vx * vx2;
        if (res > 0) {
            positive_cnt++;
        }
    }

    if (valid_cnt < 1) return NAN;
    return positive_cnt * 1.0 / valid_cnt;
}

double dummy_dcor(const std::vector<double>& x_, const std::vector<double>& y_) {
    std::vector<double> x, y;
    for (size_t j = 1; j < x_.size(); ++j) {
        x.push_back(x_[j] - x_[j - 1]);
        y.push_back(y_[j] - y_[j - 1]);
    }
    long double axx = 0, ayy = 0, axy = 0;
    std::size_t nv = 0;
    for (size_t i = 0; i < x.size(); ++i) {
        if (std::isnan(x[i]) || std::isnan(y[i])) continue;
        nv++;
        axx += x[i] * x[i];
        ayy += y[i] * y[i];
        axy += x[i] * y[i];
    }
    if (nv < 1) {
        return NAN;
    } else {
        long double d = axx * ayy;
        return d > 0 ? axy / sqrt(d) : NAN;
    }
}

double dummy_r2(const std::vector<double>& x_, const std::vector<double>& y_, double a, double b, int window) {
    double err = 0, diff = 0;
    double mean_y = mean(y_);
    int valid_count = 0;
    for (int i = 0; i < window; ++i) {
        if (std::isfinite(x_[i]) && std::isfinite(y_[i])) {
            ++valid_count;
            double fit = a + b * x_[i];
            double res = y_[i] - fit;
            err += res * res;
            diff += (y_[i] - mean_y) * (y_[i] - mean_y);
        }
    }
    if (valid_count > 1)
        return 1. - err / diff;
    else
        return NAN;
}

double dummy_r2(const std::vector<double>& x1_, const std::vector<double>& x2_, const std::vector<double>& y_,
                double b0, double b1, double b2, int window) {
    double err = 0, diff = 0;
    double mean_y = mean(y_);
    int valid_count = 0;
    for (int i = 0; i < window; ++i) {
        if (std::isfinite(x1_[i]) && std::isfinite(x2_[i]) && std::isfinite(y_[i])) {
            ++valid_count;
            double fit = b0 + b1 * x1_[i] + b2 * x2_[i];
            double res = y_[i] - fit;
            err += res * res;
            diff += (y_[i] - mean_y) * (y_[i] - mean_y);
        }
    }
    if (valid_count > 1)
        return 1. - err / diff;
    else
        return NAN;
}

double dummy_r2_no_slope(const std::vector<double>& x_, const std::vector<double>& y_, double b, int window) {
    double err = 0, diff = 0;
    int valid_count = 0;
    for (int i = 0; i < window; ++i) {
        if (std::isfinite(x_[i]) && std::isfinite(y_[i])) {
            ++valid_count;
            double fit = b * x_[i];
            double res = y_[i] - fit;
            err += res * res;
            diff += y_[i] * y_[i];
        }
    }
    if (valid_count > 0)
        return 1. - err / diff;
    else
        return NAN;
}

double dummy_r2_no_slope(const std::vector<double>& x1_, const std::vector<double>& x2_, const std::vector<double>& y_,
                         double b1, double b2, int window) {
    double err = 0, diff = 0;
    int valid_count = 0;
    for (int i = 0; i < window; ++i) {
        if (std::isfinite(x1_[i]) && std::isfinite(x2_[i]) && std::isfinite(y_[i])) {
            ++valid_count;
            double fit = b1 * x1_[i] + b2 * x2_[i];
            double res = y_[i] - fit;
            err += res * res;
            diff += y_[i] * y_[i];
        }
    }
    if (valid_count > 1)
        return 1. - err / diff;
    else
        return NAN;
}

double dummy_skew(const vector<double>& data_) {
    double mean_ = mean(data_);
    double std_ = 0;
    int valid_count = 0;
    for (double i : data_) {
        if (std::isfinite(i)) {
            ++valid_count;
            std_ += std::pow((i - mean_), 2);
        }
    }
    if (valid_count < 2) return NAN;
    std_ = sqrt(std_ / valid_count);
    if (std_ < 1e-7) return NAN;
    double ret = 0;
    for (double i : data_) {
        if (std::isfinite(i)) {
            ret += std::pow((i - mean_) / std_, 3);
        }
    }
    return ret / valid_count;
}

double dummy_kurtosis(const std::vector<double>& data_) {
    double mean_ = mean(data_);
    double std_ = 0;
    int valid_count = 0;
    for (double i : data_) {
        if (std::isfinite(i)) {
            ++valid_count;
            std_ += std::pow((i - mean_), 2);
        }
    }
    if (valid_count < 2) return NAN;
    std_ = sqrt(std_ / valid_count);
    if (std_ < 1e-7) return NAN;
    double ret = 0;
    for (double i : data_) {
        if (std::isfinite(i)) {
            ret += std::pow((i - mean_) / std_, 4);
        }
    }
    return ret / valid_count - 3.0;
}
}  // namespace ornate
