#include <math_dummy.h>
#include <math_vector.h>
#include <algorithm>
#include <cmath>
#include <cstddef>

namespace ornate {
int dummy_add(int x, int y) { return x + y; }
int universal_answer() { return 42; }
double dummy_ts_cross(const double *x_, const double *y, std::size_t i, std::size_t n) {
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

double dummy_ts_argmax(const double *x, std::size_t i, std::size_t n) {
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
    if (max_idx < 0) return NAN;
    return N > 1 ? (1 - max_idx / (N - 1.0)) : 0;
}

double dummy_ts_backward_cpn(const double *x, std::size_t i, std::size_t n, int sign) {
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

double dummy_ts_acp(const double *x, int i, int n, int lag) {
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

double dummy_dcor(const std::vector<double> &x_, const std::vector<double> &y_) {
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
    if (nv < 3) {
        return NAN;
    } else {
        long double d = axx * ayy;
        return d > 0 ? axy / sqrt(d) : NAN;
    }
}

double dummy_r2(const std::vector<double> &x_, const std::vector<double> &y_, double a, double b, int window) {
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

double dummy_r2(const std::vector<double> &x1_, const std::vector<double> &x2_, const std::vector<double> &y_,
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

double dummy_r2_no_slope(const std::vector<double> &x_, const std::vector<double> &y_, double b, int window) {
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

double dummy_r2_no_slope(const std::vector<double> &x1_, const std::vector<double> &x2_, const std::vector<double> &y_,
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

namespace dminner {
double mean(const double *x, int start, int end, double fill) {
    double res = 0, n = 0;
    for (int i = start; i < end; ++i) {
        if (std::isfinite(x[i])) {
            res += x[i];
            n++;
        }
    }
    return n > 0 ? res / n : fill;
}

double quantile(const double *x, double q, int start, int end, double fill) {
    double res;
    std::vector<double> y;
    for (int i = start; i < end; ++i) {
        if (std::isfinite(x[i])) {
            y.push_back(x[i]);
        }
    }

    std::size_t ny = y.size();
    if (ny == 0) {
        return fill;
    }

    double idx = (ny - 1) * q;
    double idx_lb = std::floor(idx);
    double idx_ub = std::ceil(idx);
    if (idx_lb == idx_ub) {
        std::nth_element(y.begin(), y.begin() + idx, y.end());
        res = y[idx];
    } else {
        std::nth_element(y.begin(), y.begin() + idx_ub, y.end());
        std::nth_element(y.begin(), y.begin() + idx_lb, y.begin() + idx_ub);
        res = y[idx_lb] * (idx_ub - idx) + y[idx_ub] * (idx - idx_lb);
    }
    return res;
}

double max(const double *x, int start, int end, double fill) {
    double res = std::numeric_limits<double>::min();
    int n = 0;
    for (int i = start; i < end; ++i) {
        if (std::isfinite(x[i]) && x[i] > res) {
            res = x[i];
            ++n;
        }
    }
    return n > 0 ? res : fill;
}

double min(const double *x, int start, int end, double fill) {
    double res = std::numeric_limits<double>::max();
    int n = 0;
    for (int i = start; i < end; ++i) {
        if (std::isfinite(x[i]) && x[i] < res) {
            res = x[i];
            ++n;
        }
    }
    return n > 0 ? res : fill;
}
}  // namespace dminner

double cond_mean(const COMPARE &g, int method, const double *x, const double *y, double q, int start, int end,
                 double fill) {
    double xx = x[end - 1], res = 0;
    if (method == 1) {
        xx = dminner::mean(x, start, end, fill);
    } else if (method == 2) {
        xx = dminner::quantile(x, q, start, end, fill);
    } else if (method == 3) {
        xx = dminner::max(x, start, end, fill);
    } else if (method == 4) {
        xx = dminner::min(x, start, end, fill);
    }
    int n = 0;
    for (int i = start; i < end; ++i) {
        if (g(x[i], xx) && std::isfinite(y[i])) {
            res += y[i];
            ++n;
        }
    }
    res = n > 0 ? res / n : fill;
    return res;
}

double cond_max(const COMPARE &g, int method, const double *x, const double *y, double q, int start, int end,
                double fill) {
    double xx = x[end - 1], res = std::numeric_limits<double>::min();
    if (method == 1) {
        xx = dminner::mean(x, start, end, fill);
    } else if (method == 2) {
        xx = dminner::quantile(x, q, start, end, fill);
    }
    for (int i = start; i < end; ++i) {
        if (g(x[i], xx) && y[i] > res) {
            res = y[i];
        }
    }
    if (res == std::numeric_limits<double>::min()) res = fill;
    return res;
}

double cond_min(const COMPARE &g, int method, const double *x, const double *y, double q, int start, int end,
                double fill) {
    double xx = x[end - 1], res = std::numeric_limits<double>::max();
    if (method == 1) {
        xx = dminner::mean(x, start, end, fill);
    } else if (method == 2) {
        xx = dminner::quantile(x, q, start, end, fill);
    }
    for (int i = start; i < end; ++i) {
        if (g(x[i], xx) && y[i] < res) {
            res = y[i];
        }
    }
    if (res == std::numeric_limits<double>::max()) res = fill;
    return res;
}

double cond_sd(const COMPARE &g, int method, const double *x, const double *y, double q, int start, int end,
               double fill) {
    double xx = x[end - 1], sx = 0, sxx = 0;
    if (method == 1) {
        xx = dminner::mean(x, start, end, fill);
    } else if (method == 2) {
        xx = dminner::quantile(x, q, start, end, fill);
    }
    int n = 0;
    for (int i = start; i < end; ++i) {
        if (g(x[i], xx) && std::isfinite(y[i])) {
            sx += y[i];
            sxx += y[i] * y[i];
            ++n;
        }
    }
    return n >= 2 ? sqrt((sxx / n - sx * sx / (n * n)) * n / (n - 1)) : fill;
}

double cond_sum(const COMPARE &g, int method, const double *x, const double *y, double q, int start, int end,
                double fill) {
    double xx = x[end - 1], res = 0;
    if (method == 1) {
        xx = dminner::mean(x, start, end, fill);
    } else if (method == 2) {
        xx = dminner::quantile(x, q, start, end, fill);
    } else if (method == 3) {
        xx = dminner::max(x, start, end, fill);
    } else if (method == 4) {
        xx = dminner::min(x, start, end, fill);
    }
    for (int i = start; i < end; ++i) {
        if (g(x[i], xx) && std::isfinite(y[i])) {
            res += y[i];
        }
    }
    return res;
}

std::vector<double> ts_cond_stat(COND_FUNC f, COMPARE g, int method, const std::vector<double> &cond,
                                 const std::vector<double> &value, int n, double q, double fill = 0, int least = 3,
                                 bool partial = true) {
    int size = cond.size();
    std::vector<double> output(size, fill);
    if (least < 3) throw std::range_error("'least' must >= 3");

    const double *x = cond.data();
    const double *y = value.data();
    double *pout = output.data();

    for (int i = (partial ? least : n) - 1; i < size; i++) {
        pout[i] = f(g, method, x, y, q, std::max(i - n + 1, 0), i + 1, fill);
    }
    return output;
}

static bool gte(double x, double y) { return x >= y; }

static bool lte(double x, double y) { return x <= y; }

std::vector<double> ts_gte_mean(const std::vector<double> &cond, const std::vector<double> &value, int n, double q,
                                double fill, int method, int least, bool partial) {
    std::vector<double> output = ts_cond_stat(cond_mean, gte, method, cond, value, n, q, fill, least, partial);
    return output;
}

std::vector<double> ts_lte_mean(const std::vector<double> &cond, const std::vector<double> &value, int n, double q,
                                double fill, int method, int least, bool partial) {
    std::vector<double> output = ts_cond_stat(cond_mean, lte, method, cond, value, n, q, fill, least, partial);
    return output;
}

std::vector<double> ts_gte_max(const std::vector<double> &cond, const std::vector<double> &value, int n, double q,
                               double fill, int method, int least, bool partial) {
    std::vector<double> output = ts_cond_stat(cond_max, gte, method, cond, value, n, q, fill, least, partial);
    return output;
}

std::vector<double> ts_lte_max(const std::vector<double> &cond, const std::vector<double> &value, int n, double q,
                               double fill, int method, int least, bool partial) {
    std::vector<double> output = ts_cond_stat(cond_max, lte, method, cond, value, n, q, fill, least, partial);
    return output;
}

std::vector<double> ts_gte_min(const std::vector<double> &cond, const std::vector<double> &value, int n, double q,
                               double fill, int method, int least, bool partial) {
    std::vector<double> output = ts_cond_stat(cond_min, gte, method, cond, value, n, q, fill, least, partial);
    return output;
}

std::vector<double> ts_lte_min(const std::vector<double> &cond, const std::vector<double> &value, int n, double q,
                               double fill, int method, int least, bool partial) {
    std::vector<double> output = ts_cond_stat(cond_min, lte, method, cond, value, n, q, fill, least, partial);
    return output;
}

std::vector<double> ts_gte_sd(const std::vector<double> &cond, const std::vector<double> &value, int n, double q,
                              double fill, int method, int least, bool partial) {
    std::vector<double> output = ts_cond_stat(cond_sd, gte, method, cond, value, n, q, fill, least, partial);
    return output;
}

std::vector<double> ts_lte_sd(const std::vector<double> &cond, const std::vector<double> &value, int n, double q,
                              double fill, int method, int least, bool partial) {
    std::vector<double> output = ts_cond_stat(cond_sd, lte, method, cond, value, n, q, fill, least, partial);
    return output;
}

std::vector<double> ts_gte_sum(const std::vector<double> &cond, const std::vector<double> &value, int n, double q,
                               double fill, int method, int least, bool partial) {
    std::vector<double> output = ts_cond_stat(cond_sum, gte, method, cond, value, n, q, fill, least, partial);
    return output;
}

std::vector<double> ts_lte_sum(const std::vector<double> &cond, const std::vector<double> &value, int n, double q,
                               double fill, int method, int least, bool partial) {
    std::vector<double> output = ts_cond_stat(cond_sum, lte, method, cond, value, n, q, fill, least, partial);
    return output;
}
}  // namespace ornate
