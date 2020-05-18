#ifndef ORNATE_MATH_STATS_ROLLING_RB_H
#define ORNATE_MATH_STATS_ROLLING_RB_H

#include "math_utils.h"

namespace ornate {

template <typename T = double>
struct rolling_rb_base {
    int window_size;
    int m_count{0}, m_valid_count{0};
    int m_head_index{0};
    std::vector<T> m_container;

    rolling_rb_base(int size) : window_size{size + 1} { m_container.resize(window_size, T()); }
};

struct rolling_mean_rb : public rolling_rb_base<double> {
    double total_sum{0};
    double mean{0};

    rolling_mean_rb(int size) : rolling_rb_base<double>(size) {}

    void delete_old() {
        int old_index = m_head_index - window_size;
        if (old_index < 0) old_index += window_size;
        const double& old_value = m_container[old_index];
        if (std::isfinite(old_value)) {
            total_sum -= old_value;
            --m_valid_count;
        }
    }
    void add_new() {
        const double& new_value = m_container[m_head_index - 1];
        if (std::isfinite(new_value)) {
            total_sum += new_value;
            ++m_valid_count;
        }

        if (m_valid_count > 0)
            mean = total_sum / m_valid_count;
        else
            mean = NAN;
    }
    double operator()(double data) {
        m_container[m_head_index++] = data;
        ++m_count;

        if (m_count >= window_size) {
            delete_old();
        }
        add_new();
        if (m_head_index == window_size) m_head_index = 0;
        return mean;
    }
};

struct rolling_variance_rb : public rolling_rb_base<double> {
    double total_sum{0}, total_square_sum{0};
    double mean{NAN}, variance{NAN};

    rolling_variance_rb(int size) : rolling_rb_base<double>(size) {}

    void delete_old() {
        int old_index = m_head_index - window_size;
        if (old_index < 0) old_index += window_size;
        const double& old_value = m_container[old_index];
        if (std::isfinite(old_value)) {
            total_sum -= old_value;
            total_square_sum -= old_value * old_value;
            --m_valid_count;
        }
    }
    void add_new() {
        const double& new_value = m_container[m_head_index - 1];
        if (std::isfinite(new_value)) {
            total_sum += new_value;
            total_square_sum += new_value * new_value;
            ++m_valid_count;
        }

        if (m_valid_count > 1) {
            mean = total_sum / m_valid_count;
            variance = (total_square_sum - mean * mean * m_valid_count) / (m_valid_count - 1);
        } else {
            mean = NAN;
            variance = NAN;
        }
    }
    double operator()(double data) {
        m_container[m_head_index++] = data;
        ++m_count;

        if (m_count >= window_size) {
            delete_old();
        }
        add_new();
        if (m_head_index == window_size) m_head_index = 0;
        return variance;
    }
};

struct rolling_cov_rb {
    double sumx{0}, sumy{0}, sumxy{0};
    double covariance{NAN};
    int window_size;
    int m_count{0}, m_valid_count{0};
    int m_head_index{0};
    std::vector<double> m_container0;
    std::vector<double> m_container1;

    rolling_cov_rb(int size) : window_size{size + 1} {
        m_container0.resize(window_size, 0);
        m_container1.resize(window_size, 0);
    }

    void delete_old() {
        int old_index = m_head_index - window_size;
        if (old_index < 0) old_index += window_size;
        double old_value0 = m_container0[old_index];
        double old_value1 = m_container1[old_index];

        if (std::isfinite(old_value0) && std::isfinite(old_value1)) {
            sumxy -= old_value0 * old_value1;
            sumx -= old_value0;
            sumy -= old_value1;
            --m_valid_count;
        }
    }
    void add_new() {
        double data0 = m_container0[m_head_index - 1];
        double data1 = m_container1[m_head_index - 1];

        if (std::isfinite(data0) && std::isfinite(data1)) {
            sumxy += data0 * data1;
            sumx += data0;
            sumy += data1;
            ++m_valid_count;
        }

        if (m_valid_count > 1) {
            double mean_x = sumx / m_valid_count;
            double mean_y = sumy / m_valid_count;
            covariance = (sumxy - mean_x * mean_y * m_valid_count) / (m_valid_count - 1);
        } else {
            covariance = NAN;
        }
    }
    double operator()(double data0, double data1) {
        m_container0[m_head_index] = data0;
        m_container1[m_head_index++] = data1;
        ++m_count;

        if (m_count >= window_size) {
            delete_old();
        }
        add_new();
        if (m_head_index == window_size) m_head_index = 0;
        return covariance;
    }
};

struct rolling_corr_rb {
    double sumx{0}, sumy{0}, sumxy{0}, sum_x2{0}, sum_y2{0};
    double corr{NAN};
    int window_size;
    int m_count{0}, m_valid_count{0};
    int m_head_index{0};
    std::vector<double> m_container0;
    std::vector<double> m_container1;

    rolling_corr_rb(int size) : window_size{size + 1} {
        m_container0.resize(window_size, 0);
        m_container1.resize(window_size, 0);
    }

    void delete_old() {
        int old_index = m_head_index - window_size;
        if (old_index < 0) old_index += window_size;
        double old_value0 = m_container0[old_index];
        double old_value1 = m_container1[old_index];

        if (std::isfinite(old_value0) && std::isfinite(old_value1)) {
            sumxy -= old_value0 * old_value1;
            sumx -= old_value0;
            sumy -= old_value1;
            sum_x2 -= old_value0 * old_value0;
            sum_y2 -= old_value1 * old_value1;
            --m_valid_count;
        }
    }
    void add_new() {
        double data0 = m_container0[m_head_index - 1];
        double data1 = m_container1[m_head_index - 1];

        if (std::isfinite(data0) && std::isfinite(data1)) {
            sumxy += data0 * data1;
            sumx += data0;
            sumy += data1;
            sum_x2 += data0 * data0;
            sum_y2 += data1 * data1;
            ++m_valid_count;
        }

        if (m_valid_count > 1) {
            double mean_x = sumx / m_valid_count;
            double mean_y = sumy / m_valid_count;
            double cov = (sumxy - mean_x * mean_y * m_valid_count) / (m_valid_count - 1);
            double var1 = (sum_x2 - mean_x * mean_x * m_valid_count) / (m_valid_count - 1);
            double var2 = (sum_y2 - mean_y * mean_y * m_valid_count) / (m_valid_count - 1);
            if (std::isfinite(var1) && std::isfinite(var2)) {
                double numerator = sqrt(var1) * sqrt(var2);
                if (numerator < epsilon)
                    corr = NAN;
                else
                    corr = cov / numerator;
            }

        } else {
            corr = NAN;
        }
    }
    double operator()(double data0, double data1) {
        m_container0[m_head_index] = data0;
        m_container1[m_head_index++] = data1;
        ++m_count;

        if (m_count >= window_size) {
            delete_old();
        }
        add_new();
        if (m_head_index == window_size) m_head_index = 0;
        return corr;
    }
};

}  // namespace ornate

#endif
