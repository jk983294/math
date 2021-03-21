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
template <typename T, typename T1>
double mean_weighted(const T *data, const T1 *weight, int32_t n) {
    double sum = 0, wgt_sum = 0;
    uint32_t count = 0;
    for (int32_t i = 0; i < n; i++) {
        if (isvalid(data[i]) && isvalid(weight[i])) {
            sum += data[i] * weight[i];
            wgt_sum += weight[i];
            count++;
        }
    }
    if (count > 0)
        return sum / wgt_sum;
    else
        return NAN;
}

/**
 * get wgt mean. will consider nan, [start_idx, end_idx), -1 to use all
 */
template <typename T, typename T1>
double mean_weighted(IN const std::vector<T> &n, IN const std::vector<T1> &n1, int32_t start_idx = -1,
                     int32_t end_idx = -1) {
    if (start_idx < 0) start_idx = 0;
    if (end_idx < 0) end_idx = static_cast<int32_t>(n.size());
    return mean_weighted(n.data() + start_idx, n1.data() + start_idx, end_idx - start_idx);
}

template <typename T, typename T1>
double sum_weighted(const T *data, const T1 *weight, int32_t n) {
    double sum = 0;
    uint32_t count = 0;
    for (int32_t i = 0; i < n; i++) {
        if (isvalid(data[i]) && isvalid(weight[i])) {
            sum += data[i] * weight[i];
            count++;
        }
    }
    if (count > 0)
        return sum;
    else
        return NAN;
}

template <typename T, typename T1>
double sum_weighted(const std::vector<T> &n, const std::vector<T1> &n1, int32_t start_idx = -1, int32_t end_idx = -1) {
    if (start_idx < 0) start_idx = 0;
    if (end_idx < 0) end_idx = static_cast<int32_t>(n.size());
    return sum_weighted(n.data() + start_idx, n1.data() + start_idx, end_idx - start_idx);
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
    std_x = std::sqrt((sum_x2 - mean_x * mean_x * count) / (count - 1));
    std_y = std::sqrt((sum_y2 - mean_y * mean_y * count) / (count - 1));
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
int __weighted_cov(const T *x, const T *y, const T *weight, size_t num, double &cov_, double &std_x, double &std_y) {
    double sum_x = 0, sum_x2 = 0, sum_xy = 0, sum_y = 0, sum_y2 = 0, sum_weight = 0;
    int count = 0;
    for (size_t i = 0; i < num; ++i) {
        if (!std::isfinite(x[i]) || !std::isfinite(y[i]) || !std::isfinite(weight[i])) continue;
        ++count;
        sum_weight += weight[i];
        sum_x += x[i] * weight[i];
        sum_y += y[i] * weight[i];
        sum_x2 += x[i] * x[i] * weight[i];
        sum_y2 += y[i] * y[i] * weight[i];
        sum_xy += x[i] * y[i] * weight[i];
    }
    if (count < 2) return count;
    double mean_x = sum_x / sum_weight;
    double mean_y = sum_y / sum_weight;
    cov_ = (sum_xy / sum_weight - mean_x * mean_y);
    std_x = std::sqrt(sum_x2 / sum_weight - mean_x * mean_x);
    std_y = std::sqrt(sum_y2 / sum_weight - mean_y * mean_y);
    return count;
}

template <typename T = float>
double weighted_corr(const T *x, const T *y, const T *weight, size_t num) {
    double cov_, std_x, std_y;
    if (weight == nullptr) {
        if (__cov(x, y, num, cov_, std_x, std_y) < 2) return NAN;
    } else {
        if (__weighted_cov(x, y, weight, num, cov_, std_x, std_y) < 2) return NAN;
    }
    if (std_x < epsilon || std_y < epsilon) return NAN;
    return cov_ / std_x / std_y;
}

template <typename T = float>
double weighted_corr(const std::vector<T> &x, const std::vector<T> &y, const std::vector<T> &weight) {
    return weighted_corr(&x[0], &y[0], &weight[0], x.size());
}

template <typename T = float>
int __rcov(const T *x, const T *y, size_t num, double &cov_, double &std_x, double &std_y, int sign) {
    double sum_x2 = 0, sum_xy = 0, sum_y2 = 0;
    int count = 0;
    for (size_t i = 0; i < num; ++i) {
        if (!std::isfinite(x[i]) || !std::isfinite(y[i])) continue;
        if (sign != 0 && y[i] * sign < 0) continue;
        ++count;
        sum_x2 += x[i] * x[i];
        sum_y2 += y[i] * y[i];
        sum_xy += x[i] * y[i];
    }
    if (count < 2) return count;
    cov_ = sum_xy;
    std_x = std::sqrt(sum_x2);
    std_y = std::sqrt(sum_y2);
    return count;
}

template <typename T = float>
double rcor(const T *x, const T *y, int num, int sign = 0) {
    double cov_, std_x, std_y;
    if (__rcov(x, y, num, cov_, std_x, std_y, sign) < 2) return NAN;
    if (std_x < epsilon || std_y < epsilon) return NAN;
    return cov_ / std_x / std_y;
}

template <typename T = float>
double rcor(const std::vector<T> &x, const std::vector<T> &y, int sign = 0) {
    return rcor(&x[0], &y[0], x.size(), sign);
}

template <typename T = float>
int __weighted_rcov(const T *x, const T *y, const T *weight, int num, double &cov_, double &std_x, double &std_y,
                    int sign) {
    double sum_x2 = 0, sum_xy = 0, sum_y2 = 0;
    int count = 0;
    for (int i = 0; i < num; ++i) {
        if (!std::isfinite(x[i]) || !std::isfinite(y[i]) || !std::isfinite(weight[i]) || weight[i] <= 0) continue;
        if (sign != 0 && y[i] * sign < 0) continue;
        ++count;
        sum_x2 += x[i] * x[i] * weight[i];
        sum_y2 += y[i] * y[i] * weight[i];
        sum_xy += x[i] * y[i] * weight[i];
    }
    if (count < 2) return count;
    cov_ = sum_xy;
    std_x = std::sqrt(sum_x2);
    std_y = std::sqrt(sum_y2);
    return count;
}

template <typename T = float>
double weighted_rcor(const T *x, const T *y, const T *weight, int num, int sign = 0) {
    double cov_, std_x, std_y;
    if (weight == nullptr) {
        if (__rcov(x, y, num, cov_, std_x, std_y, sign) < 2) return NAN;
    } else {
        if (__weighted_rcov(x, y, weight, num, cov_, std_x, std_y, sign) < 2) return NAN;
    }
    if (std_x < epsilon || std_y < epsilon) return NAN;
    return cov_ / std_x / std_y;
}

template <typename T = float>
double weighted_rcor(const std::vector<T> &x, const std::vector<T> &y, const std::vector<T> &weight, int sign = 0) {
    return weighted_rcor(&x[0], &y[0], &weight[0], x.size(), sign);
}

template <typename T1 = float, typename T2 = float>
int __sign_cov(const T1 *x, const T2 *y, size_t num, double &cov_, double &std_x, double &std_y) {
    double sum_x = 0, sum_x2 = 0, sum_xy = 0, sum_y = 0, sum_y2 = 0;
    int count = 0;
    for (size_t i = 0; i < num; ++i) {
        double v1 = x[i];
        double v2 = y[i];
        if (!std::isfinite(v1) || !std::isfinite(v2)) continue;

        if (v1 > 0) {
            v1 = 1;
            if (v2 >= 0)
                v2 = 1;
            else
                v2 = -1;
        } else if (v1 < 0) {
            v1 = -1.;
            if (v2 <= 0)
                v2 = -1;
            else
                v2 = 1;
        } else {
            if (v2 > 0) {
                v1 = v2 = 1;
            } else if (v2 < 0) {
                v1 = v2 = -1;
            } else
                continue;
        }
        ++count;
        sum_x += v1;
        sum_y += v2;
        sum_x2 += v1 * v1;
        sum_y2 += v2 * v2;
        sum_xy += v1 * v2;
    }
    if (count < 2) return count;
    double mean_x = sum_x / count;
    double mean_y = sum_y / count;
    cov_ = (sum_xy - mean_x * mean_y * count) / (count - 1);
    std_x = std::sqrt((sum_x2 - mean_x * mean_x * count) / (count - 1));
    std_y = std::sqrt((sum_y2 - mean_y * mean_y * count) / (count - 1));
    return count;
}

template <typename T1 = float, typename T2 = float>
double sign_corr(const T1 *x, const T2 *y, size_t num) {
    double cov_, std_x, std_y;
    if (__sign_cov(x, y, num, cov_, std_x, std_y) < 2) return NAN;
    if (std_x < epsilon || std_y < epsilon) return NAN;
    return cov_ / std_x / std_y;
}

template <typename T1 = float, typename T2 = float>
double sign_corr(const std::vector<T1> &x, const std::vector<T2> &y) {
    return sign_corr(&x[0], &y[0], x.size());
}

template <typename T = float>
double acf(const T *x, int num, int shift) {
    double cov_, std_x, std_y;
    if (__cov(x, x + shift, num - shift, cov_, std_x, std_y) < 2) return NAN;
    if (std_x < epsilon || std_y < epsilon) return NAN;
    return cov_ / std_x / std_y;
}

template <typename T = float>
double acf(const std::vector<T> &x, int shift) {
    return acf(x.data(), (int)x.size(), shift);
}

template <typename T = float>
int acf_half_life(const T *x, int max_life, int num, int shift_bulk = 1) {
    int left = 0, right = max_life;
    double right_acf = acf(x, num, shift_bulk * max_life);
    if (std::isfinite(right_acf) && right_acf >= 0.5) return max_life;
    for (;;) {
        int middle = (left + right) / 2;
        if (middle == left || middle == right) return middle;
        double middle_acf = acf(x, num, shift_bulk * middle);
        if (std::isfinite(right_acf)) {
            if (middle_acf > 0.5)
                left = middle;
            else
                right = middle;
        } else {
            right = middle;
        }
    }
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
    return std::sqrt((double)sum / (count - 1));
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
    if (data.empty()) return NAN;
    if (data.size() < (size_t)n) {
        int new_n = (int)data.size();
        int diff = n - i;
        if (new_n - diff >= 0) {
            return ema_decay(data.data(), new_n, new_n - diff, decay);
        } else
            return NAN;
    } else {
        return ema_decay(data.data(), n, i, decay);
    }
}

template <typename T = float>
double ema_hl(const T *data, int n, int i, double hl) {
    double decay = ema_hl2decay(hl);
    return ema_decay(data, n, i, decay);
}

template <typename T = float>
double ema_hl(const std::vector<T> &data, int n, int i, double hl) {
    double decay = ema_hl2decay(hl);
    return ema_decay(data, n, i, decay);
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

template <typename T>
double math_skew(const T *data_, int num) {
    double mean_ = mean(data_, num);
    double std_ = 0;
    int valid_count = 0;
    for (int i = 0; i < num; ++i) {
        if (std::isfinite(data_[i])) {
            ++valid_count;
            std_ += std::pow((data_[i] - mean_), 2);
        }
    }
    if (valid_count < 2) return NAN;
    std_ = sqrt(std_ / valid_count);
    if (std_ < 1e-7) return NAN;
    double ret = 0;
    for (int i = 0; i < num; ++i) {
        if (std::isfinite(data_[i])) {
            ret += std::pow((data_[i] - mean_) / std_, 3);
        }
    }
    return ret / valid_count;
}

template <typename T>
double math_kurtosis(const T *data_, int num) {
    double mean_ = mean(data_, num);
    double std_ = 0;
    int valid_count = 0;
    for (int i = 0; i < num; ++i) {
        if (std::isfinite(data_[i])) {
            ++valid_count;
            std_ += std::pow((data_[i] - mean_), 2);
        }
    }
    if (valid_count < 2) return NAN;
    std_ = sqrt(std_ / valid_count);
    if (std_ < 1e-7) return NAN;
    double ret = 0;
    for (int i = 0; i < num; ++i) {
        if (std::isfinite(data_[i])) {
            ret += std::pow((data_[i] - mean_) / std_, 4);
        }
    }
    return ret / valid_count - 3.0;
}

template <typename T>
double math_skew(const vector<T> &data_) {
    return math_skew(data_.data(), data_.size());
}

template <typename T>
double math_kurtosis(const std::vector<T> &data_) {
    return math_kurtosis(data_.data(), data_.size());
}

template <typename T = double>
int calc_na_count(const T *x, int num) {
    int cnt = 0;
    for (int i = 0; i < num; ++i) {
        if (!std::isfinite(x[i])) cnt++;
    }
    return cnt;
}

template <typename T, typename T1>
double calc_r_square(const T *y, const T1 *y_hat, int num) {
    double mean_y = 0;
    int valid_num = 0;
    for (int i = 0; i < num; i++) {
        if (!isfinite(y[i]) || !isfinite(y_hat[i])) continue;
        mean_y += y[i];
        valid_num++;
    }

    if (valid_num > 3) {
        mean_y = mean_y / valid_num;
        double ssreg = 0, sstot = 0;
        for (int i = 0; i < num; i++) {
            if (!isfinite(y[i]) || !isfinite(y_hat[i])) continue;
            ssreg += std::pow((y[i] - y_hat[i]), 2);
            sstot += std::pow((y[i] - mean_y), 2);
        }
        return 1. - ssreg / sstot;
    } else
        return NAN;
}

template <typename T>
std::pair<int, int> math_sign_count(const T *data_, int num) {
    int neg_cnt = 0, pos_cnt = 0;
    for (int i = 0; i < num; ++i) {
        if (isvalid(data_[i])) {
            if (data_[i] > 0.)
                ++pos_cnt;
            else if (data_[i] < 0.)
                ++neg_cnt;
        }
    }
    return {neg_cnt, pos_cnt};
}

template <typename T>
std::pair<double, double> math_sign_ratio(const T *data_, int num) {
    int neg_cnt, pos_cnt;
    std::tie(neg_cnt, pos_cnt) = math_sign_count(data_, num);

    if (neg_cnt + pos_cnt > 0) {
        double total_cnt = neg_cnt + pos_cnt;
        return {neg_cnt / total_cnt, pos_cnt / total_cnt};
    }
    return {NAN, NAN};
}

template <typename T>
std::pair<double, double> math_sign_ratio(const vector<T> &data_) {
    return math_sign_ratio(data_.data(), data_.size());
}

template <typename T>
double math_majority_sign_ratio(const T *data_, int num) {
    double r1, r2;
    std::tie(r1, r2) = math_sign_ratio(data_, num);
    if (std::isfinite(r1)) return std::max(r1, r2);
    return NAN;
}

template <typename T>
double math_majority_sign_ratio(const vector<T> &data_) {
    return math_majority_sign_ratio(data_.data(), data_.size());
}

template <typename T>
double cap_within_sd(const T x, double mean_, double sd, double range = 3.0) {
    if (!isvalid(x))
        return mean_;
    else if (x < mean_ - range * sd)
        return mean_ - range * sd;
    else if (x > mean_ + range * sd)
        return mean_ + range * sd;
    return x;
}

struct skip_corr_stat {
    double sum_x{0}, sum_y{0};
    double sum_x2{0}, sum_y2{0}, sum_xy{0};
    int count{0};

    skip_corr_stat() = default;

    void operator()(double x, double y) {
        if (std::isfinite(x) && std::isfinite(y)) {
            sum_x += x;
            sum_y += y;
            sum_x2 += x * x;
            sum_y2 += y * y;
            sum_xy += x * y;
            ++count;
        }
    }

    double get_corr() const {
        if (count < 2) return NAN;
        double mean_x = sum_x / count;
        double mean_y = sum_y / count;
        double cov_ = (sum_xy - mean_x * mean_y * count) / (count - 1);
        double std_x = std::sqrt((sum_x2 - mean_x * mean_x * count) / (count - 1));
        double std_y = std::sqrt((sum_y2 - mean_y * mean_y * count) / (count - 1));
        if (std_x < epsilon || std_y < epsilon) return NAN;
        return cov_ / std_x / std_y;
    }
};

template <typename T = float>
std::vector<double> skip_corr(const T *x, const T *y, int num, int skip) {
    vector<double> ret(skip, NAN);
    vector<skip_corr_stat> stats(skip);
    int idx = 0;
    for (int i = 0; i < num; ++i) {
        (stats[idx++])(x[i], y[i]);
        if (idx == skip) idx = 0;
    }
    for (int i = 0; i < skip; ++i) {
        ret[i] = stats[i].get_corr();
    }
    return ret;
}

template <typename T = float>
std::vector<double> skip_corr(const std::vector<T> &x, const std::vector<T> &y, int skip) {
    return skip_corr(x.data(), y.data(), x.size(), skip);
}

template <typename T>
vector<vector<double>> calc_histogram_stats(vector<double> cuts, const vector<T> &y_hat, const vector<T> &y) {
    vector<vector<double>> ret;
    int cut_cnt = (int)cuts.size();
    cuts.push_back(NAN);
    vector<vector<double>> histogram_y(cut_cnt + 1);
    vector<vector<double>> histogram_y_hat(cut_cnt + 1);
    ret.resize(cut_cnt + 1);
    for (size_t i = 0; i < y_hat.size(); ++i) {
        if (!isvalid(y_hat[i]) or !isvalid(y[i])) continue;

        int idx = std::lower_bound(cuts.begin(), cuts.end(), y_hat[i]) - cuts.begin() + 1;
        if (idx == 1 && y_hat[i] < cuts.front())
            idx = 0;
        else if (idx == cut_cnt + 1)
            idx = cut_cnt;
        histogram_y[idx].push_back(y[i]);
        histogram_y_hat[idx].push_back(y_hat[i]);
    }

    for (int i = 0; i < cut_cnt + 1; ++i) {
        ret[i].push_back(cuts[i]);
        ret[i].push_back(histogram_y[i].size());
        ret[i].push_back(ornate::mean(histogram_y[i]));
        ret[i].push_back(ornate::mean(histogram_y_hat[i]));
        ret[i].push_back(ornate::corr(histogram_y[i], histogram_y_hat[i]));
    }
    return ret;
}

template <typename T>
vector<vector<double>> split_histogram(vector<double> cuts, const vector<T> &y) {
    vector<vector<double>> histograms;
    int cut_cnt = (int)cuts.size();
    histograms.resize(cut_cnt + 1);
    for (size_t i = 0; i < y.size(); ++i) {
        if (!isvalid(y[i])) continue;

        int idx = std::lower_bound(cuts.begin(), cuts.end(), y[i]) - cuts.begin() + 1;
        if (idx == 1 && y[i] < cuts.front())
            idx = 0;
        else if (idx == cut_cnt + 1)
            idx = cut_cnt;
        histograms[idx].push_back(y[i]);
    }
    cuts.push_back(NAN);
    return histograms;
}

template <typename TGroup, typename T>
vector<vector<double>> split_histogram_by(vector<TGroup> cuts, const vector<TGroup> &group, const vector<T> &y) {
    vector<vector<double>> histograms;
    int cut_cnt = (int)cuts.size();
    histograms.resize(cut_cnt + 1);
    for (size_t i = 0; i < y.size(); ++i) {
        if (!isvalid(y[i])) continue;

        int idx = std::lower_bound(cuts.begin(), cuts.end(), group[i]) - cuts.begin() + 1;
        if (idx == 1 && group[i] < cuts.front())
            idx = 0;
        else if (idx == cut_cnt + 1)
            idx = cut_cnt;
        histograms[idx].push_back(y[i]);
    }
    cuts.push_back(get_nan<TGroup>());
    return histograms;
}

template <typename T>
vector<int> calc_histogram_cnt(vector<double> cuts, const vector<T> &y) {
    vector<vector<double>> histograms = split_histogram(cuts, y);
    vector<int> cnts(histograms.size());
    for (const auto &s : histograms) {
        cnts.push_back((int)s.size());
    }
    return cnts;
}

template <typename T>
T median(const T *x, std::vector<T> &data, int num) {
    int valid_cnt = 0;
    for (int i = 0; i < num; ++i) {
        if (isvalid(x[i])) {
            data.push_back(x[i]);
            ++valid_cnt;
        }
    }
    if (valid_cnt > 0) {
        int n = valid_cnt / 2;
        nth_element(data.begin(), data.begin() + n, data.end());
        if (valid_cnt % 2 == 0) {
            auto max_it = std::max_element(data.begin(), data.begin() + n);
            return (*max_it + data[n]) / 2.0;
        } else {
            return data[n];
        }
    } else {
        return get_nan<T>();
    }
}

template <typename T>
T median(const T *x, int num) {
    std::vector<T> data;
    data.reserve(num);
    return median(x, data, num);
}

template <typename T>
T median(const std::vector<T> &x) {
    return median(x.data(), x.size());
}
}  // namespace ornate

#endif
