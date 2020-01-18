#ifndef ORNATE_MATH_STATS_H
#define ORNATE_MATH_STATS_H

#include "math_utils.h"

using std::isfinite;

/**
 * dof: degree of freedom, 1 means sample statistics
 * R: return type
 * T: data type
 */

namespace ornate {

template <typename T>
float mean(IN const std::vector<T> &n, int32_t start_idx = -1, int32_t end_idx = -1) {
    T ret = 0;
    if (start_idx < 0) start_idx = 0;
    if (end_idx < 0) end_idx = static_cast<int32_t>(n.size());
    uint32_t count = 0;
    for (int32_t i = start_idx; i < end_idx; ++i) {
        if (isvalid(n[i])) {
            ret += n[i];
            count++;
        }
    }
    if (count > 0)
        return ret / (float)count;
    else
        return NAN;
}

template <typename T>
T mean(const T *data, size_t n) {
    T sum = 0;
    size_t count = 0;
    for (size_t i = 0; i < n; ++i) {
        if (IsValidData(data[i])) {
            ++count;
            sum += data[i];
        }
    }
    if (count == 0)
        return NAN;
    else
        return sum / count;
}

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

//! welford's algorithm
template <int dof = 1, typename R = double, typename T>
R variance_online(R &variance, R &mean, const T data, size_t &n) {
    R delta, SSE;
    if (n == 0)
        return NAN;  // throw exception, should > 0
    else if (n > 1)
        SSE = variance * (n - 1 - dof);
    else if (n == 1) {
        variance = 0;
        if (!std::isfinite(data) || std::isnan(data)) return variance;
        mean = data;
        ++n;
        return variance;
    }

    if (std::isfinite(data) && !std::isnan(data)) {
        delta = data - mean;
        mean = mean + delta / n;
        SSE = SSE + delta * (data - mean);
        variance = SSE / (n - dof);
        ++n;
        return variance;
    } else {
        return variance;
    }
}

// adapted version of welford's algorithm
template <int dof = 1, typename R = double, typename T = double>
struct variance_rolling {
    R SSE{0}, mean{0}, variance{0};
    long count{0};
    std::deque<T> m_data;
    size_t window_size;

    explicit variance_rolling(size_t size) : window_size{size} {}

    R operator()(T data) {
        if (!std::isfinite(data) || std::isnan(data)) {  // invalid data
            m_data.push_back(data);
            if (m_data.size() > window_size) {
                if (std::isfinite(m_data[0]) && !std::isnan(m_data[0])) {
                    mean = (mean * count - m_data[0]) / (count - 1);
                    SSE -= (1.0 * (count - 1) / count) * (m_data[0] - mean) * (m_data[0] - mean);
                    m_data.pop_front();
                    --count;
                } else {
                    m_data.pop_front();
                }
            }
            if (count < 2) variance = NAN;
            variance = SSE / (count - dof);
            return variance;
        }

        if (count == 0) {  // first data
            mean = data;
            SSE = 0;
            m_data.push_back(data);
            ++count;
            if (m_data.size() > window_size)  // if the deque was full of NAN by window_size before push_back
                m_data.pop_front();
            variance = NAN;
            return variance;
        }

        m_data.push_back(data);
        ++count;

        if (m_data.size() <= window_size) {  // window not full
            SSE += (1.0 * (count - 1) / (count)) * (data - mean) * (data - mean);
            mean = mean + 1.0 * (data - mean) / count;
        } else {  // window full, now starts to roll
            if (std::isfinite(m_data[0]) && !std::isnan(m_data[0])) {
                auto old = m_data[0];
                auto old_mean = mean;
                mean += 1.0 * (data - old) / (count - 1);
                SSE +=
                    (1.0 * (count - 1) / count) * ((data - old_mean) * (data - old_mean) - (old - mean) * (old - mean));
                m_data.pop_front();
                --count;
            } else {
                // if oldest data is NAN, this is not rolling either
                SSE += (1.0 * (count - 1) / count) * (data - mean) * (data - mean);
                mean = mean + 1.0 * (data - mean) / (count);
                m_data.pop_front();
            }
        }
        if (count < 2) variance = NAN;
        variance = SSE / (count - dof);
        return variance;
    }
};

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

template <int dof = 1, typename R = double, typename T1, typename T2>
R covariance_online(R &covariance, R &meanA, R &meanB, const T1 dataA, const T2 dataB, size_t &n) {
    R deltaA, deltaB;
    R CM;
    if (n == 0) return NAN;  // throw ("cannot set n to 0");
    if (n > 1)
        CM = covariance * (n - 1 - dof);
    else if (n == 1) {
        covariance = 0;
        if (!std::isfinite(dataA) || std::isnan(dataA) || !std::isfinite(dataB) || std::isnan(dataB)) return covariance;
        meanA = dataA;
        meanB = dataB;
        ++n;
        return covariance;
    }
    if (std::isfinite(dataA) && !std::isnan(dataA) && std::isfinite(dataB) && !std::isnan(dataB)) {
        deltaA = dataA - meanA;
        meanA = meanA + deltaA / n;
        deltaB = dataB - meanB;
        meanB = meanB + deltaB / n;
        CM = CM + (n - 1) * 1.0 / n * (deltaA) * (deltaB);
        covariance = CM / (n - dof);
        ++n;
        return covariance;
    } else {
        return covariance;  // do nothing is data is not finite or nan
    }
}

// adapted version of welford's algorithm
template <int dof = 1, typename R = double, typename TA = double, typename TB = double>
struct covariance_rolling {
    R CM{0};
    R covariance{0};
    R meanA{0}, meanB{0};
    long count{0};  // only count for valid data
    std::deque<R> m_dataA;
    std::deque<R> m_dataB;
    size_t window_size;

    explicit covariance_rolling(size_t size) : window_size{size} {}

    R operator()(TA dataA, TB dataB) {
        if (!IsValidData(dataA) || !IsValidData(dataB)) {
            // NAN data also pushed in, but not count
            m_dataA.push_back(dataA);
            m_dataB.push_back(dataB);
            if (m_dataA.size() > window_size) {
                if (IsValidData(m_dataA[0]) && IsValidData(m_dataB[0])) {
                    meanA = (meanA * count - m_dataA[0]) / (count - 1);
                    meanB = (meanB * count - m_dataB[0]) / (count - 1);
                    CM -= (1.0 * (count - 1) / count) * (m_dataA[0] - meanA) * (m_dataB[0] - meanB);
                    m_dataA.pop_front();
                    m_dataB.pop_front();
                    --count;
                } else {
                    m_dataA.pop_front();
                    m_dataB.pop_front();
                }
            }
            if (count < 2)
                covariance = NAN;
            else
                covariance = CM / (count - dof);
            return covariance;
        }

        if (count == 0) {
            meanA = dataA;
            meanB = dataB;
            CM = 0;
            m_dataA.push_back(dataA);
            m_dataB.push_back(dataB);
            ++count;

            if (m_dataA.size() > window_size) {  // if the deque was full of NAN by window_size before push_back
                m_dataA.pop_front();
                m_dataB.pop_front();
            }
            covariance = NAN;
            return covariance;
        }

        m_dataA.push_back(dataA);
        m_dataB.push_back(dataB);
        ++count;

        if (m_dataA.size() <= window_size) {  // window not full
            CM += (1.0 * (count - 1) / (count)) * (dataA - meanA) * (dataB - meanB);
            meanA = meanA + 1.0 * (dataA - meanA) / count;
            meanB = meanB + 1.0 * (dataB - meanB) / count;
        } else {  // window full, now starts to roll
            if (IsValidData(m_dataA[0]) && IsValidData(m_dataB[0])) {
                auto oldA = m_dataA[0];
                auto oldB = m_dataB[0];
                auto oldmeanA = meanA;
                auto oldmeanB = meanB;
                meanA += 1.0 * (dataA - oldA) / (count - 1);
                meanB += 1.0 * (dataB - oldB) / (count - 1);
                CM += (1.0 * (count - 1) / count) *
                      ((dataA - oldmeanA) * (dataB - oldmeanB) - (oldA - meanA) * (oldB - meanB));
                m_dataA.pop_front();
                m_dataB.pop_front();
                --count;
            } else {
                // if oldest data is NAN, this is not rolling either
                CM += (1.0 * (count - 1) / count) * (dataA - meanA) * (dataB - meanB);
                meanA = meanA + 1.0 * (dataA - meanA) / (count);
                meanB = meanB + 1.0 * (dataB - meanB) / (count);
                m_dataA.pop_front();
                m_dataB.pop_front();
            }
        }
        if (count < 2)
            covariance = NAN;
        else
            covariance = CM / (count - dof);
        return covariance;
    }
};

template <typename T = float>
size_t __cov(const T *x, const T *y, size_t num, T &cov_, T &std_x, T &std_y) {
    T sum_x = 0, sum_x2 = 0, sum_xy = 0, sum_y = 0, sum_y2 = 0;
    size_t count = 0;
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
    T mean_x = sum_x / count;
    T mean_y = sum_y / count;
    cov_ = sum_xy / count - mean_x * mean_y;
    std_x = sqrtf(sum_x2 / count - mean_x * mean_x);
    std_y = sqrtf(sum_y2 / count - mean_y * mean_y);
    return count;
}

struct OnlineCorrelation {
    double m_correlation{0};
    double meanA{0}, meanB{0};
    double meanAA{0}, meanBB{0};
    size_t n{1}, nA{1}, nB{1};
    double varA{0}, varB{0}, cov{0};

    void reset() {
        m_correlation = 0;
        meanA = 0;
        meanB = 0;
        meanAA = 0;
        meanBB = 0;
        n = 1;
        nA = 1;
        nB = 1;
        varA = 0;
        varB = 0;
        cov = 0;
    }

    void Push(double dataA, double dataB) {
        if (!isfinite(dataA) || !isfinite(dataB)) {
            dataA = dataB = NAN;
        }
        ornate::variance_online(varA, meanA, dataA, nA);
        ornate::variance_online(varB, meanB, dataB, nB);
        ornate::covariance_online(cov, meanAA, meanBB, dataA, dataB, n);
    }
    double Result() {
        double numerator = sqrt(varA) * sqrt(varB);
        if (numerator < 1e-7)
            m_correlation = NAN;
        else
            m_correlation = cov / numerator;
        return m_correlation;
    }
};

template <typename T = float>
T corr(const T *x, const T *y, size_t num) {
    T cov_, std_x, std_y;
    if (__cov(x, y, num, cov_, std_x, std_y) < 2) return NAN;
    if (std_x < 0.000001 || std_y < 0.000001) return NAN;
    return cov_ / std_x / std_y;
}

template <typename T = float>
T corr(const std::vector<T> &x, const std::vector<T> &y) {
    return corr(&x[0], &y[0], x.size());
}

template <typename T = float>
T cov(const T *x, const T *y, size_t num) {
    T cov_, std_x, std_y;
    if (__cov(x, y, num, cov_, std_x, std_y) < 2) return NAN;
    if (std_x < 0.000001 || std_y < 0.000001) return NAN;
    return cov_;
}

template <typename T = float>
T cov(const std::vector<T> &x, const std::vector<T> &y) {
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
