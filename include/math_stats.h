#ifndef ORNATE_MATH_STATS_H
#define ORNATE_MATH_STATS_H

#include "math_stats_rolling.h"
#include "math_utils.h"
#include "math_vector.h"

using std::isfinite;

/**
 * dof: degree of freedom, 1 means sample statistics
 * R: return type
 * T: data type
 */

namespace ornate {

template <typename R = double, typename T1, typename T2>
R mean_weighted(const T1 *data, const T2 *weight, size_t n) {
    R sum = 0;
    size_t count = 0;
    for (size_t i = 0; i < n; ++i) {
        if (IsValidData(data[i])) {
            ++count;
            sum += data[i] * weight[i];
        }
    }
    if (count == 0)
        return NAN;
    else
        return sum / count;
}

template <int dof = 1, typename R = double, typename T>
R variance_two_pass(const T *data, size_t n) {
    R sum1 = 0, sum2 = 0, sum3 = 0, mean = 0;
    for (size_t i = 0; i < n; ++i) {
        sum1 = sum1 + data[i];
    }
    mean = sum1 / n;

    for (size_t i = 0; i < n; ++i) {
        sum2 = sum2 + (data[i] - mean) * (data[i] - mean);
        sum3 = sum3 + (data[i] - mean);
    }
    return (sum2 - sum3 * sum3 / n) / (n - dof);
}

template <int dof = 1, typename R = double, typename T>
R variance_one_pass(const T *data, size_t n) {
    R mean = 0, M2 = 0, delta;
    size_t count = 0;
    for (size_t i = 0; i < n; ++i) {
        if (std::isfinite(data[i])) {
            ++count;
            delta = data[i] - mean;
            mean = mean + delta / count;
            M2 = M2 + delta * (data[i] - mean);
        }
    }
    if (count == 0)
        return NAN;
    else if (count < 2)
        return 0;
    return M2 / (count - dof);
}

template <int dof = 1, typename R = double, typename T, typename T2>
R variance_weighted(const T *data, const T2 *weight, size_t n) {
    T2 sum_weight = 0, temp, ratio;
    R mean = 0, M2 = 0, delta;
    size_t count = 0;
    for (size_t i = 0; i < n; ++i) {
        if (std::isfinite(data[i])) {
            ++count;
            temp = weight[i] + sum_weight;
            delta = data[i] - mean;
            ratio = delta * weight[i] / temp;
            mean = mean + ratio;
            M2 = M2 + delta * ratio * sum_weight;
            sum_weight = temp;
        }
    }
    if (count == 0)
        return NAN;
    else if (count < 2)
        return 0;
    else
        return M2 / sum_weight * n / (n - dof);
}

template <int dof = 1, typename R = double, typename T>
R variance(const T *data, size_t n) {
    return variance_one_pass<dof, R, T>(data, n);
}

template <int dof = 1, typename R = double, typename T1, typename T2>
R covariance_two_pass(const T1 *data1, const T2 *data2, size_t n) {
    R sum1 = 0, sum2 = 0, mean1 = 0, mean2 = 0, covariance = 0;
    size_t count = 0;
    for (size_t i = 0; i < n; ++i) {
        if (std::isfinite(data1[i]) && std::isfinite(data2[i])) {
            ++count;
            sum1 = sum1 + data1[i];
            sum2 = sum2 + data2[i];
        }
    }
    mean1 = sum1 / count;
    mean2 = sum2 / count;
    for (size_t i = 0; i < n; ++i) {
        if (std::isfinite(data1[i]) && std::isfinite(data1[i])) {
            covariance += (data1[i] - mean1) * (data2[i] - mean2) / (count - dof);
        }
    }
    return covariance;
}

template <int dof = 1, typename R = double, typename T1, typename T2>
R covariance_one_pass(const T1 *data1, const T2 *data2, size_t n) {
    R mean1 = 0, mean2 = 0, C12 = 0, delta1 = 0, delta2 = 0;
    size_t count = 0;
    for (size_t i = 0; i < n; ++i, ++count) {
        if (std::isfinite(data1[i]) && std::isfinite(data2[i])) {
            delta1 = (data1[i] - mean1) / (count + 1);
            delta2 = (data2[i] - mean2) / (count + 1);
            C12 += (data1[i] - mean1) * (data2[i] - mean2) * count / (count + 1);
            mean1 += delta1;
            mean2 += delta2;
        }
    }
    if (count == 0)
        return NAN;
    else if (count < 2)
        return 0;
    return C12 / (count - dof);
}

template <int dof = 1, typename R = double, typename T1, typename T2>
R covariance(const T1 *data1, const T2 *data2, size_t n) {
    return covariance_one_pass<dof, R, T1, T2>(data1, data2, n);
}

template <typename T = float>
int __cov(const T *x, const T *y, size_t num, double &cov_, double &std_x, double &std_y) {
    double sum_x = 0, sum_x2 = 0, sum_xy = 0, sum_y = 0, sum_y2 = 0;
    int count = 0;
    for (size_t i = 0; i < num; ++i) {
        if (!std::isfinite(x[i]) || !std::isfinite(y[i])) continue;
        ++count;
        sum_x += x[i];
        sum_y += y[i];
        sum_x2 += x[i] * x[i];
        sum_y2 += y[i] * y[i];
        sum_xy += x[i] * y[i];
    }
    if (count < 2) return count;
    double mean_x = sum_x / count;
    double mean_y = sum_y / count;
    cov_ = (sum_xy - mean_x * mean_y * count) / (count - 1);
    std_x = sqrtf((sum_x2 - mean_x * mean_x * count) / (count - 1));
    std_y = sqrtf((sum_y2 - mean_y * mean_y * count) / (count - 1));
    return count;
}

template <typename T = float>
double corr(const T *x, const T *y, size_t num) {
    double cov_, std_x, std_y;
    if (__cov(x, y, num, cov_, std_x, std_y) < 2) return NAN;
    if (std_x < epsilon || std_y < epsilon) return NAN;
    return cov_ / std_x / std_y;
}

template <typename T = float>
double corr(const std::vector<T> &x, const std::vector<T> &y) {
    return corr(&x[0], &y[0], x.size());
}

template <typename T = float>
double cov(const T *x, const T *y, size_t num) {
    double cov_, std_x, std_y;
    if (__cov(x, y, num, cov_, std_x, std_y) < 2) return NAN;
    if (std_x < epsilon || std_y < epsilon) return NAN;
    return cov_;
}

template <typename T = float>
double cov(const std::vector<T> &x, const std::vector<T> &y) {
    return cov(&x[0], &y[0], x.size());
}

template <typename T = float>
T std(const T *data, int32_t n) {
    T sum = 0;
    int32_t count = 0;
    for (int32_t i = 0; i < n; i++) {
        if (std::isfinite(data[i])) {
            sum += data[i];
            ++count;
        }
    }
    if (count == 0) return NAN;
    T mean_ = sum / count;
    sum = 0;
    for (int32_t i = 0; i < n; i++) {
        if (std::isfinite(data[i])) {
            sum += (data[i] - mean_) * ((data[i] - mean_));
        }
    }
    return sqrtf(sum / (count - 1));
}

template <typename T = float>
T std(const std::vector<T> &data, int32_t start_idx = -1, int32_t end_idx = -1) {
    if (start_idx < 0) start_idx = 0;
    if (end_idx < 0) end_idx = static_cast<int32_t>(data.size());
    return std(&data[start_idx], end_idx - start_idx);
}

template <typename T = float>
bool _regression(const T *y, const T *x, size_t num, OUT T *a, OUT T *b, OUT T *R) {
    std::vector<T> ny(num, 0);
    std::vector<T> nx(num, 0);
    size_t valid_num = 0;
    for (size_t i = 0; i < num; i++) {
        if (!isfinite(y[i]) || !isfinite(x[i])) continue;
        ny[valid_num] = y[i];
        nx[valid_num] = x[i];
        valid_num++;
    }
    if (valid_num < 2) {
        *a = NAN;
        *b = NAN;
        return false;
    }
    T mx = mean(nx.data(), valid_num);
    T my = mean(ny.data(), valid_num);
    T denominator = 0, numerator = 0;
    for (size_t i = 0; i < valid_num; i++) denominator += (nx[i] - mx) * (nx[i] - mx);
    for (size_t i = 0; i < valid_num; i++) numerator += (ny[i] - my) * (nx[i] - mx);
    *b = numerator / denominator;
    *a = my - *b * mx;
    if (R != nullptr && valid_num > 3) {
        T err = 0;
        for (size_t i = 0; i < valid_num; i++) {
            T t = ny[i] - (*a + *b * nx[i]);
            err += t * t;
        }
        T dis = 0;
        for (size_t i = 0; i < valid_num; i++) {
            dis += (ny[i] - my) * (ny[i] - my);
        }
        *R = 1 - (err / (valid_num - 2)) / (dis / (valid_num - 1));
    }
    return true;
}

/**
 * find the regression of two series. y = a + b * x, calc a, b, and R(if not nullptr)
 * return false when data is not sufficient
 */
template <typename T = float>
bool regression(const std::vector<T> &y, const std::vector<T> &x, OUT T *a, OUT T *b, OUT T *R = nullptr) {
    return _regression(&y[0], &x[0], x.size(), a, b, R);
}

/**
 * y = b0 + b1 * x1 + b2 * x2
 */
template <typename T = float>
inline bool regression3(const std::vector<T> &y, const std::vector<T> &x1, const std::vector<T> &x2, OUT T *b0,
                        OUT T *b1, OUT T *b2) {
    size_t valid_num = 0;
    double sum_x12 = 0, sum_x1_2 = 0, sum_x2_2 = 0, sum_y_2 = 0, sum_x1 = 0, sum_x2 = 0, sum_y = 0, sum_x1y = 0,
           sum_x2y = 0;
    for (size_t i = 0; i < y.size(); i++) {
        if (!isfinite(y[i]) || !isfinite(x1[i]) || !isfinite(x2[i])) continue;
        sum_x12 += x1[i] * x2[i];
        sum_x1_2 += x1[i] * x1[i];
        sum_x2_2 += x2[i] * x2[i];
        sum_y_2 += y[i] * y[i];
        sum_x1y += x1[i] * y[i];
        sum_x2y += x2[i] * y[i];
        sum_x1 += x1[i];
        sum_x2 += x2[i];
        sum_y += y[i];
        valid_num++;
    }

    if (valid_num < 3) {
        *b0 = NAN;
        *b1 = NAN;
        *b2 = NAN;
        return false;
    }
    double sum_X1_2 = sum_x1_2 - sum_x1 * sum_x1 / valid_num;
    double sum_X2_2 = sum_x2_2 - sum_x2 * sum_x2 / valid_num;
    double sum_X1Y = sum_x1y - sum_x1 * sum_y / valid_num;
    double sum_X2Y = sum_x2y - sum_x2 * sum_y / valid_num;
    double sum_X12 = sum_x12 - sum_x1 * sum_x2 / valid_num;

    double denominator = sum_X1_2 * sum_X2_2 - sum_X12 * sum_X12;
    double b1_ = (sum_X2_2 * sum_X1Y - sum_X12 * sum_X2Y) / denominator;
    double b2_ = (sum_X1_2 * sum_X2Y - sum_X12 * sum_X1Y) / denominator;
    double b0_ = (sum_y - b1_ * sum_x1 - b2_ * sum_x2) / valid_num;
    *b0 = b0_;
    *b1 = b1_;
    *b2 = b2_;
    return true;
}

template <typename T = float>
T _ols(const T *y, const T *x, size_t num) {
    T sum_x2 = 0, sum_xy = 0;
    size_t valid_num = 0;
    for (size_t i = 0; i < num; i++) {
        if (!isfinite(y[i]) || !isfinite(x[i])) continue;
        sum_xy += y[i] * x[i];
        sum_x2 += x[i] * x[i];
        valid_num++;
    }
    if (valid_num < 1) {
        return NAN;
    }
    return sum_xy / sum_x2;
}
/**
 * ols differ from regression is that it has no intercept, i.e y = b * x
 */
template <typename T = float>
T ols(const std::vector<T> &y, const std::vector<T> &x) {
    return _ols(&y[0], &x[0], x.size());
}

/**
 * ols differ from regression is that it has no intercept, i.e y = b1 * x1 + b2 * x2
 */
template <typename T = float>
inline bool ols(const std::vector<T> &y, const std::vector<T> &x1, const std::vector<T> &x2, OUT T *b1, OUT T *b2) {
    size_t valid_num = 0;
    long double sum_x12 = 0, sum_x1_2 = 0, sum_x2_2 = 0, sum_x1y = 0, sum_x2y = 0;
    for (size_t i = 0; i < y.size(); i++) {
        if (!isfinite(y[i]) || !isfinite(x1[i]) || !isfinite(x2[i])) continue;
        sum_x12 += x1[i] * x2[i];
        sum_x1_2 += x1[i] * x1[i];
        sum_x2_2 += x2[i] * x2[i];
        sum_x1y += x1[i] * y[i];
        sum_x2y += x2[i] * y[i];
        valid_num++;
    }

    if (valid_num < 2) {
        *b1 = NAN;
        *b2 = NAN;
        return false;
    }
    long double denominator = sum_x12 * sum_x12 - sum_x1_2 * sum_x2_2;
    long double b1_ = (sum_x2y * sum_x12 - sum_x1y * sum_x2_2) / denominator;
    long double b2_ = (sum_x1y * sum_x12 - sum_x2y * sum_x1_2) / denominator;
    *b1 = b1_;
    *b2 = b2_;
    return true;
}

template <typename T = float>
T quantile(T *data, int size, double q) {
    auto itr = std::partition(data, data + size, [](auto i) { return std::isfinite(i); });
    int valid_count = itr - data;
    if (valid_count <= 0) return NAN;

    double idx = (valid_count - 1) * q;
    long nth_lb = std::lround(std::floor(idx));
    long nth_ub = std::lround(std::ceil(idx));
    if (nth_lb < 0) nth_lb = 0;
    if (nth_lb >= valid_count) nth_lb = valid_count - 1;
    if (nth_ub < 0) nth_ub = 0;
    if (nth_ub >= valid_count) nth_ub = valid_count - 1;
    if (nth_lb == nth_ub) {
        std::nth_element(data, data + nth_lb, data + valid_count);
        return data[nth_lb];
    } else {
        std::nth_element(data, data + nth_ub, data + valid_count);
        std::nth_element(data, data + nth_lb, data + nth_ub);
        return data[nth_lb] * (nth_ub - idx) + data[nth_ub] * (idx - nth_lb);
    }
}

template <typename T = float>
T quantile(std::vector<T> &data, double q) {
    return quantile(data.data(), data.size(), q);
}

template <typename T = float>
double normal_ema(const std::vector<T> &data, int num) {
    double alpha = 2.0f / static_cast<double>(1 + num);
    double ret = data.front();
    for (size_t i = 1; i < data.size(); ++i) {
        ret = ret * (1 - alpha) + data[i] * alpha;
    }
    return ret;
}

template <typename T = float>
double ema_decay(const T *data, int n, int i, double decay) {
    if (i < n - 1) return NAN;
    double v = 0, w = 0;
    for (int ii = 0; ii < n; ++ii) {
        if (std::isfinite(data[i - ii])) {
            v += data[i - ii] * pow(decay, ii);
            w += pow(decay, ii);
        }
    }
    return v / w;
}

template <typename T = float>
double ema_decay(const std::vector<T> &data, int n, int i, double decay) {
    if ((size_t)n > data.size()) return NAN;
    return ema_decay(data.data(), n, i, decay);
}

template <typename T = float>
double ema_hl(const T *data, int n, int i, double hl) {
    double decay = ema_hl2decay(hl);
    return ema_decay(data, n, i, decay);
}

template <typename T = float>
double ema_hl(const std::vector<T> &data, int n, int i, double hl) {
    if ((size_t)n > data.size()) return NAN;
    return ema_hl(data.data(), n, i, hl);
}

template <typename T = float>
T _slope_no_intercept(const T *y, size_t num) {
    T sum_x2 = 0, sum_xy = 0;
    size_t valid_num = 0;
    for (size_t i = 0; i < num; i++) {
        if (!isfinite(y[i])) continue;
        sum_xy += y[i] * i;
        sum_x2 += i * i;
        valid_num++;
    }
    if (valid_num < 1) {
        return NAN;
    }
    return sum_xy / sum_x2;
}

template <typename T = float>
T slope_no_intercept(const std::vector<T> &y) {
    return _slope_no_intercept(&y[0], y.size());
}

template <typename T = float>
T _ts_slope(const T *y, size_t num) {
    T sum_x2 = 0, sum_xy = 0, sum_x{0}, sum_y{0};
    size_t valid_num = 0;
    for (size_t i = 0; i < num; i++) {
        if (!isfinite(y[i])) continue;
        sum_xy += y[i] * i;
        sum_x2 += i * i;
        sum_x += i;
        sum_y += y[i];
        valid_num++;
    }
    if (valid_num > 0) {
        long double cov_xy = sum_xy * valid_num - sum_x * sum_y;
        long double var_x = sum_x2 * valid_num - sum_x * sum_x;
        return var_x > 0 ? cov_xy / var_x : NAN;
    }
    return NAN;
}

template <typename T = float>
T ts_slope(const std::vector<T> &y) {
    return _ts_slope(&y[0], y.size());
}

template <typename T = float>
T _ts_sharpe(const T *y, size_t num) {
    T sum_y2 = 0, sum_y{0};
    size_t valid_num = 0;
    for (size_t i = 0; i < num; i++) {
        if (!isfinite(y[i])) continue;
        sum_y2 += y[i] * y[i];
        sum_y += y[i];
        valid_num++;
    }

    if (valid_num <= 1) return NAN;
    long double mean = sum_y / valid_num;
    long double sd = sqrt((sum_y2 - mean * mean * valid_num) / (valid_num - 1.0));
    return sd > 0 ? mean / sd : NAN;
}

template <typename T = float>
T ts_sharpe(const std::vector<T> &y) {
    return _ts_sharpe(&y[0], y.size());
}

template <typename T = float>
T _ts_scale(const T *y, size_t num) {
    T sum_abs_y{0};
    size_t valid_num = 0;
    for (size_t i = 0; i < num; i++) {
        if (!isfinite(y[i])) continue;
        sum_abs_y += std::abs(y[i]);
        valid_num++;
    }

    T new_value = y[num - 1];
    if (std::isfinite(new_value)) {
        double abs_mean = sum_abs_y / valid_num;
        return abs_mean > 1e-6 ? (new_value / abs_mean) : NAN;
    } else {
        return NAN;
    }
}

template <typename T = float>
T ts_scale(const std::vector<T> &y) {
    return _ts_scale(&y[0], y.size());
}
}  // namespace ornate

#endif
