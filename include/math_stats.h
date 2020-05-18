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
}  // namespace ornate

#endif
