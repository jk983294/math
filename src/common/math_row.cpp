#include <math_row.h>
#include <algorithm>
#include <functional>
#include <limits>
#include <utility>
#include <vector>

namespace ornate {

double row_sum(const double* a, int len) {
    double total_sum{0};
    int valid_count{0};
    for (int i = 0; i < len; ++i) {
        auto val = a[i];
        if (std::isfinite(val)) {
            total_sum += val;
            ++valid_count;
        }
    }
    if (valid_count > 0)
        return total_sum;
    else
        return NAN;
}

double row_sum2(const double* a, const double* b, int len) {
    double total_sum{0};
    int valid_count{0};
    for (int i = 0; i < len; ++i) {
        auto val0 = a[i];
        auto val1 = b[i];
        if (std::isfinite(val0) && std::isfinite(val1)) {
            total_sum += val0 * val1;
            ++valid_count;
        }
    }
    if (valid_count > 0)
        return total_sum;
    else
        return NAN;
}

int row_cut(const double* a, double threshold, int len) {
    for (int i = 0; i < len; ++i) {
        auto val = a[i];
        if (!std::isfinite(val)) {  // when nan found, cut return nan
            return -1;
        }

        if (val >= threshold) {
            return i + 1;
        } else if (i + 1 >= len) {
            return len;
        }
    }
    return -1;
}

int row_sum_cut(const double* a, double threshold, int len) {
    int idx{-1};
    double _sum{0};
    for (int i = 0; i < len; ++i) {
        auto val = a[i];
        if (!std::isfinite(val)) {  // when nan found, cut return nan
            return -1;
        }

        _sum += val;
        if (_sum >= threshold) {
            return i + 1;
        } else if (i + 1 >= len) {
            return len;
        }
    }
    if (idx > 0)
        return idx;
    else
        return -1;
}

double row_mean(const double* a, int len) {
    double total_sum{0};
    int valid_count{0};
    for (int i = 0; i < len; ++i) {
        auto val = a[i];
        if (std::isfinite(val)) {
            total_sum += val;
            ++valid_count;
        }
    }
    if (valid_count > 0)
        return total_sum / valid_count;
    else
        return NAN;
}

double row_mean2(const double* a, const double* b, int len) {
    double total_sum{0};
    double total_weight{0};
    int valid_count{0};
    for (int i = 0; i < len; ++i) {
        auto val0 = a[i];
        auto val1 = b[i];
        if (std::isfinite(val0) && std::isfinite(val1)) {
            total_sum += val0 * val1;
            total_weight += val1;
            ++valid_count;
        }
    }
    if (valid_count > 0)
        return total_sum / total_weight;
    else
        return NAN;
}

double row_ema(const double* a, int len, double hl) {
    const double decay{1. / hl};
    double v{0};
    int valid_count{0};
    for (int i = 0; i < len; ++i) {
        auto val0 = a[i];
        if (std::isfinite(val0)) {
            if (valid_count == 0)
                v = val0;
            else
                v += decay * (val0 - v);
            ++valid_count;
        }
    }
    if (valid_count > 0)
        return v;
    else
        return NAN;
}

double row_ema2(const double* a, const double* b, int len, double hl) {
    const double decay{1. / hl};
    double v{0};
    int valid_count{0};
    for (int i = 0; i < len; ++i) {
        auto val0 = a[i];
        auto weight_ = (i == 0 ? 0.0 : b[i]);
        if (std::isfinite(val0) && std::isfinite(weight_)) {
            if (valid_count == 0)
                v = val0;
            else
                v += decay * weight_ * (val0 - v);
            ++valid_count;
        }
    }
    if (valid_count > 0)
        return v;
    else
        return NAN;
}

double row_sd(const double* a, int len) {
    double sx{0}, sx2{0};
    int valid_count{0};
    for (int i = 0; i < len; ++i) {
        auto val = a[i];
        if (std::isfinite(val)) {
            sx += val;
            sx2 += val * val;
            ++valid_count;
        }
    }
    if (valid_count > 1)
        return std::sqrt(sx2 / (valid_count - 1) - sx * sx / (valid_count * (valid_count - 1)));
    else
        return NAN;
}

double row_sd2(const double* a, const double* b, int len) {
    double sx{0}, sx2{0};
    double sw{0};
    int valid_count{0};
    for (int i = 0; i < len; ++i) {
        auto val0 = a[i];
        auto weight_ = b[i];
        if (std::isfinite(val0) && std::isfinite(weight_)) {
            sx += val0 * weight_;
            sx2 += val0 * val0 * weight_;
            sw += weight_;
            ++valid_count;
        }
    }
    if (sw != 0 && valid_count > 1)
        return std::sqrt((sx2 / sw - sx * sx / (sw * sw)) * valid_count / (valid_count - 1));
    else
        return NAN;
}

double row_median(const double* a, int len) {
    std::vector<double> v;
    v.reserve(len);
    for (int i = 0; i < len; ++i) {
        if (std::isfinite(a[i])) {
            v.emplace_back(a[i]);
        }
    }
    if (v.empty()) return NAN;

    std::sort(v.begin(), v.end());
    if (v.size() % 2 == 0) {
        return (v[v.size() / 2 - 1] + v[v.size() / 2]) / 2.0;
    } else {
        return v[v.size() / 2];
    }
}

double row_slope(const double* a, int len, bool xnorm, bool ynorm, bool intercept) {
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    int valid_count{0};
    for (int i = 0; i < len; ++i) {
        auto val = a[i];
        if (std::isfinite(val)) {
            double x = i + 1.;
            sx += x;
            sy += val;
            sxx += x * x;
            sxy += x * val;
            ++valid_count;
        }
    }

    if (valid_count < 2) return NAN;
    if (xnorm) {
        sxx = sxx / (sx * sx);
        sxy = sxy / sx;
        sx = 1;
    }
    if (ynorm) {
        sxy = sxy / sy;
        sy = 1;
    }
    if (intercept) {
        long double varx = sxx / valid_count - sx * sx / (valid_count * valid_count);
        if (varx > 1e-12) {
            long double covxy = sxy / valid_count - sx * sy / (valid_count * valid_count);
            return covxy / varx;
        } else {
            return NAN;
        }
    } else {
        return sxx > 1e-6 ? sxy / sxx : NAN;
    }
}

double row_slope2(const double* a, const double* b, int len, bool xnorm, bool ynorm, bool intercept) {
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    int valid_count{0};
    for (int i = 0; i < len; ++i) {
        if (std::isfinite(a[i]) && std::isfinite(b[i])) {
            sx += a[i];
            sy += b[i];
            sxx += a[i] * a[i];
            sxy += a[i] * b[i];
            ++valid_count;
        }
    }

    if (valid_count < 2) return NAN;
    if (xnorm) {
        sxx = sxx / (sx * sx);
        sxy = sxy / sx;
        sx = 1;
    }
    if (ynorm) {
        sxy = sxy / sy;
        sy = 1;
    }
    if (intercept) {
        long double varx = sxx / valid_count - sx * sx / (valid_count * valid_count);
        if (varx > 1e-12) {
            long double covxy = sxy / valid_count - sx * sy / (valid_count * valid_count);
            return covxy / varx;
        } else {
            return NAN;
        }
    } else {
        return sxx > 1e-6 ? sxy / sxx : NAN;
    }
}

double row_diff_mad(const double* a, int len) {
    double _sum = 0;
    int valid_count{0};
    for (int i = 1; i < len; ++i) {
        if (std::isfinite(a[i - 1]) && std::isfinite(a[i])) {
            _sum += fabs(a[i] - a[i - 1]);
            ++valid_count;
        }
    }
    return valid_count > 0 ? _sum / valid_count : NAN;
}

double row_sum_ad(const double* a, int len, double center) {
    double _sum = 0;
    int valid_count{0};
    for (int i = 0; i < len; ++i) {
        if (std::isfinite(a[i])) {
            _sum += fabs(a[i] - center);
            ++valid_count;
        }
    }
    return valid_count > 0 ? _sum : NAN;
}

double row_sum_ad2(const double* a, const double* b, int len, double center) {
    double _sum = 0, _sum_w{0};
    int valid_count{0};
    for (int i = 0; i < len; ++i) {
        if (std::isfinite(a[i]) && std::isfinite(b[i])) {
            _sum += b[i] * fabs(a[i] - center);
            _sum_w += b[i];
            ++valid_count;
        }
    }
    return valid_count > 0 && _sum_w > 0 ? _sum / _sum_w : NAN;
}

double row_sum_ed(const double* a, int len, double center, double lambda) {
    double _sum = 0;
    int valid_count{0};
    for (int i = 0; i < len; ++i) {
        if (std::isfinite(a[i])) {
            _sum += exp(-lambda * fabs(a[i] - center));
            ++valid_count;
        }
    }
    return valid_count > 0 ? _sum : NAN;
}

double row_sum_ed2(const double* a, const double* b, int len, double center, double lambda) {
    double _sum = 0, _sum_w{0};
    int valid_count{0};
    for (int i = 0; i < len; ++i) {
        if (std::isfinite(a[i]) && std::isfinite(b[i])) {
            _sum += b[i] * exp(-lambda * fabs(a[i] - center));
            _sum_w += b[i];
            ++valid_count;
        }
    }
    return valid_count > 0 && _sum_w > 0 ? _sum / _sum_w : NAN;
}

double row_max(const double* a, int len) {
    double max_val{NAN};
    for (int i = 0; i < len; ++i) {
        auto val = a[i];
        if (std::isfinite(val)) {
            if (std::isnan(max_val) || val > max_val) {
                max_val = val;
            }
        }
    }
    return max_val;
}

double row_min(const double* a, int len) {
    double min_val{NAN};
    for (int i = 0; i < len; ++i) {
        auto val = a[i];
        if (std::isfinite(val)) {
            if (std::isnan(min_val) || val < min_val) {
                min_val = val;
            }
        }
    }
    return min_val;
}

int row_which_max(const double* a, int len) {
    double max_val{NAN};
    int idx{-1};
    for (int i = 0; i < len; ++i) {
        auto val = a[i];
        if (std::isfinite(val)) {
            if (std::isnan(max_val) || val > max_val) {
                max_val = val;
                idx = i + 1;
            }
        }
    }
    return idx;
}

int row_which_min(const double* a, int len) {
    double min_val{NAN};
    int idx{-1};
    for (int i = 0; i < len; ++i) {
        auto val = a[i];
        if (std::isfinite(val)) {
            if (std::isnan(min_val) || val < min_val) {
                min_val = val;
                idx = i + 1;
            }
        }
    }
    return idx;
}

template <typename TCmp>
int row_which_op(const double* a, double val, int len) {
    for (int i = 0; i < len; ++i) {
        auto data = a[i];
        if (!std::isfinite(data)) {  // when nan found, return nan
            return -1;
        }
        if (TCmp()(data, val)) {
            return i + 1;
        }
    }
    return -1;
}

int row_which_eq(const double* a, double val, int len) { return row_which_op<std::equal_to<double>>(a, val, len); }
int row_which_gt(const double* a, double val, int len) { return row_which_op<std::greater<double>>(a, val, len); }
int row_which_lt(const double* a, double val, int len) { return row_which_op<std::less<double>>(a, val, len); }
int row_which_gte(const double* a, double val, int len) {
    return row_which_op<std::greater_equal<double>>(a, val, len);
}
int row_which_lte(const double* a, double val, int len) { return row_which_op<std::less_equal<double>>(a, val, len); }
double row_extract(const double* a, int idx, int len) {
    if (idx <= 0 || idx > len) {
        return NAN;
    } else {
        return a[idx - 1];
    }
}

}  // namespace mathrow