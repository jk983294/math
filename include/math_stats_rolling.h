#ifndef ORNATE_MATH_STATS_ROLLING_H
#define ORNATE_MATH_STATS_ROLLING_H

#include "math_stats_rolling_rb.h"
#include "math_utils.h"
#include "math_vector.h"

using std::isfinite;

/**
 * dof: degree of freedom, 1 means sample statistics
 * R: return type
 * T: data type
 */

namespace ornate {

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
            if (count < 2)
                variance = NAN;
            else
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

struct OnlineCorrelation {
    double m_correlation{0};
    double meanA{0}, meanB{0};
    double meanAA{0}, meanBB{0};
    size_t n{1}, nA{1}, nB{1};
    double varA{0}, varB{0}, cov{0};

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

struct corr_rolling {
    size_t window_size;
    covariance_rolling<> cov_rolling;
    variance_rolling<> var_rolling_a;
    variance_rolling<> var_rolling_b;

    explicit corr_rolling(size_t size)
        : window_size{size}, cov_rolling(size), var_rolling_a(size), var_rolling_b(size) {}

    double operator()(double dataA, double dataB) {
        if (!isfinite(dataA) || !isfinite(dataB)) {
            dataA = dataB = NAN;
        }
        auto cov = cov_rolling(dataA, dataB);
        auto var1 = var_rolling_a(dataA);
        auto var2 = var_rolling_b(dataB);
        if (std::isfinite(var1) && std::isfinite(var2)) {
            double numerator = sqrt(var1) * sqrt(var2);
            if (numerator < epsilon)
                return NAN;
            else
                return cov / numerator;
        } else {
            return NAN;
        }
    }
};

template <int dof = 1, typename R = double, typename T = double>
struct skew_rolling {
    R mean{0};
    long count{0};
    std::deque<T> m_data;
    size_t window_size;

    explicit skew_rolling(size_t size) : window_size{size} {}

    R operator()(T data) {
        if (std::isfinite(data)) {  // invalid data
            m_data.push_back(data);
            ++count;

            if (m_data.size() <= window_size) {  // window not full
                mean = mean + 1.0 * (data - mean) / count;
            } else {  // window full, now starts to roll
                mean += 1.0 * (data - m_data[0]) / (count - 1);
                m_data.pop_front();
                --count;
            }
        }

        if (count < 2) return NAN;
        R a = 0, b = 0;
        for (auto d : m_data) {
            a += std::pow(d - mean, 3);
            b += (d - mean) * (d - mean);
        }
        return (a / count) / (std::pow(b / (count - dof), 1.5));
    }
};

template <typename R = double, typename T = double>
struct mean_rolling {
    R mean{0};
    long count{0};
    std::deque<T> m_data;
    size_t window_size;

    explicit mean_rolling(size_t size) : window_size{size} {}

    R operator()(T data) {
        if (std::isfinite(data)) {  // valid data
            m_data.push_back(data);
            ++count;

            if (m_data.size() <= window_size) {  // window not full
                mean = mean + 1.0 * (data - mean) / count;
            } else {  // window full, now starts to roll
                mean += 1.0 * (data - m_data[0]) / (count - 1);
                m_data.pop_front();
                --count;
            }
        }

        if (count < 1) return NAN;
        return mean;
    }
};

}  // namespace ornate

#endif
