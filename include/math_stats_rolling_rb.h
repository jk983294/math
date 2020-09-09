#ifndef ORNATE_MATH_STATS_ROLLING_RB_H
#define ORNATE_MATH_STATS_ROLLING_RB_H

#include <algorithm>
#include <functional>
#include "math_utils.h"

namespace ornate {

template <typename T = double>
struct rolling_rb_base {
    int window_size;
    int m_count{0}, m_valid_count{0};
    int m_head_index{0};
    std::vector<T> m_container;

    rolling_rb_base(int size) : window_size{size + 1} { m_container.resize(window_size, T()); }

    int get_old_index() {
        int old_index = m_head_index - window_size;
        if (old_index < 0) old_index += window_size;
        return old_index;
    }
};

struct rolling_mean_rb : public rolling_rb_base<double> {
    double total_sum{0};
    double mean{0};

    rolling_mean_rb(int size) : rolling_rb_base<double>(size) {}

    void delete_old() {
        int old_index = get_old_index();
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
        int old_index = get_old_index();
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

struct rolling_skew_rb : public rolling_rb_base<double> {
    double total_x1{0}, total_x2{0}, total_x3{0};
    double skew{0};

    rolling_skew_rb(int size) : rolling_rb_base<double>(size) {}

    void delete_old() {
        int old_index = get_old_index();
        const double& old_value = m_container[old_index];
        if (std::isfinite(old_value)) {
            total_x1 -= old_value;
            total_x2 -= old_value * old_value;
            total_x3 -= old_value * old_value * old_value;
            --m_valid_count;
        }
    }
    void add_new() {
        const double& new_value = m_container[m_head_index - 1];
        if (std::isfinite(new_value)) {
            total_x1 += new_value;
            total_x2 += new_value * new_value;
            total_x3 += new_value * new_value * new_value;
            ++m_valid_count;
        }

        if (m_valid_count >= 2) {
            double mean = total_x1 / m_valid_count;
            double mean2 = mean * mean;
            double var = total_x2 / m_valid_count - mean2;
            if (var <= 1e-14)
                skew = NAN;
            else {
                double mean3 = mean * mean * mean;
                double m3 = total_x3 / m_valid_count - 3 * mean * total_x2 / m_valid_count + 2 * mean3;
                skew = m3 / std::pow(var, 1.5);
            }
        } else
            skew = NAN;
    }
    double operator()(double data) {
        m_container[m_head_index++] = data;
        ++m_count;

        if (m_count >= window_size) {
            delete_old();
        }
        add_new();
        if (m_head_index == window_size) m_head_index = 0;
        return skew;
    }
};

struct rolling_kurtosis_rb : public rolling_rb_base<double> {
    double total_x1{0}, total_x2{0}, total_x3{0}, total_x4{0};
    double kurtosis{0};

    rolling_kurtosis_rb(int size) : rolling_rb_base<double>(size) {}

    void delete_old() {
        int old_index = get_old_index();
        const double& old_value = m_container[old_index];
        if (std::isfinite(old_value)) {
            total_x1 -= old_value;
            total_x2 -= old_value * old_value;
            total_x3 -= old_value * old_value * old_value;
            total_x4 -= old_value * old_value * old_value * old_value;
            --m_valid_count;
        }
    }
    void add_new() {
        const double& new_value = m_container[m_head_index - 1];
        if (std::isfinite(new_value)) {
            total_x1 += new_value;
            total_x2 += new_value * new_value;
            total_x3 += new_value * new_value * new_value;
            total_x4 += new_value * new_value * new_value * new_value;
            ++m_valid_count;
        }

        if (m_valid_count >= 2) {
            double mean = total_x1 / m_valid_count;
            double mean2 = mean * mean;
            double var = total_x2 / m_valid_count - mean2;
            if (var <= 1e-14)
                kurtosis = NAN;
            else {
                double mean3 = mean2 * mean;
                double mean4 = mean3 * mean;
                double m4 = total_x4 / m_valid_count - 4 * mean * total_x3 / m_valid_count +
                            6 * mean2 * total_x2 / m_valid_count - 3 * mean4;
                kurtosis = m4 / std::pow(var, 2) - 3.0;
            }
        } else
            kurtosis = NAN;
    }
    double operator()(double data) {
        m_container[m_head_index++] = data;
        ++m_count;

        if (m_count >= window_size) {
            delete_old();
        }
        add_new();
        if (m_head_index == window_size) m_head_index = 0;
        return kurtosis;
    }
};

template <typename T>
struct rolling_rank_base_rb {
    struct SortItem {
        T data = T();
        int seq{-1};
        SortItem() = default;
        explicit SortItem(T data_, int seq_) : data{data_}, seq{seq_} {}
        bool operator<(const SortItem& a) const { return data < a.data; }
        bool operator==(const SortItem& a) const { return data == a.data && seq == a.seq; }
        bool operator<(const T a) const { return data < a; }
        bool operator==(const T a) const { return data == a; }
    };

    int window_size;
    int m_count{0}, m_valid_count{0};
    int m_head_index{0};
    std::vector<T> m_container;
    std::vector<SortItem> m_sorted_data;

    explicit rolling_rank_base_rb(int size) : window_size{size + 1} {
        m_container.resize(window_size, T());
        m_sorted_data.resize(size);
    }

    int get_old_index() {
        int old_index = m_head_index - window_size;
        if (old_index < 0) old_index += window_size;
        return old_index;
    }

    /**
     * upper_bound 向左找到比data大的位置，然后 rotate 把 m_sorted_data[idx] 挪到该位置
     * 每次找比data大的位置插入，这样还保证了seq维度上的stable insert
     * @return idx where new data reside
     */
    int insert_left(T data, int idx) {
        m_sorted_data[idx].data = data;
        m_sorted_data[idx].seq = m_count - 1;
        auto itr = std::upper_bound(m_sorted_data.begin(), m_sorted_data.begin() + idx, m_sorted_data[idx]);
        std::rotate(itr, m_sorted_data.begin() + idx, m_sorted_data.begin() + idx + 1);
        return itr - m_sorted_data.begin();
    }

    /**
     * upper_bound 向右找到比data大的位置，然后 rotate 把 m_sorted_data[idx] 挪到该位置
     * 每次找比data大的位置插入，这样还保证了seq维度上的stable insert
     * @return idx where new data reside
     */
    int insert_right(T data, int idx) {
        m_sorted_data[idx].data = data;
        m_sorted_data[idx].seq = m_count - 1;
        auto itr =
            std::upper_bound(m_sorted_data.begin() + idx, m_sorted_data.begin() + m_valid_count, m_sorted_data[idx]);
        std::rotate(m_sorted_data.begin() + idx, m_sorted_data.begin() + idx + 1, itr);
        return int(itr - m_sorted_data.begin()) - 1;
    }

    int find_insert_point(T data, int seq) {
        auto itr = std::lower_bound(m_sorted_data.begin(), m_sorted_data.begin() + m_valid_count, data);
        int pos = itr - m_sorted_data.begin();
        while (pos < m_valid_count && m_sorted_data[pos].seq != seq && m_sorted_data[pos].data == data) pos++;
        return pos;
    }

    /**
     * 因为新插入的值是在相等值的最右边(stable sort)，向前找到地一个不等于它的位置 + 1就是 lower bound
     */
    int lower_bound(T data, int idx) {
        while (idx >= 0 && m_sorted_data[idx].data == data) --idx;
        return idx + 1;
    }

    std::pair<int, int> handle(T new_value) {
        m_container[m_head_index++] = new_value;
        ++m_count;

        int idx = -1, lower = -1;
        if (m_count >= window_size) {
            int old_index = get_old_index();
            T old_value = m_container[old_index];
            if (isvalid(old_value)) {
                int insert_pos_index = find_insert_point(old_value, m_count - window_size);

                if (isvalid(new_value)) {
                    if (new_value > old_value) {
                        idx = insert_right(new_value, insert_pos_index);
                    } else {
                        idx = insert_left(new_value, insert_pos_index);
                    }

                    lower = lower_bound(new_value, idx);
                } else {
                    // as new_value is NAN, move back data one step forward
                    std::rotate(m_sorted_data.begin() + insert_pos_index, m_sorted_data.begin() + insert_pos_index + 1,
                                m_sorted_data.begin() + m_valid_count);
                    --m_valid_count;
                }
            } else {
                if (isvalid(new_value)) {
                    idx = insert_left(new_value, m_valid_count);
                    lower = lower_bound(new_value, idx);
                    m_valid_count++;
                }
            }
        } else {  // 数据还不够，直接 insert sort
            if (isvalid(new_value)) {
                idx = insert_left(new_value, m_valid_count);
                lower = lower_bound(new_value, idx);
                m_valid_count++;
            }
        }

        if (m_head_index == window_size) m_head_index = 0;
        return {lower, idx};
    }
};

template <typename T>
struct rolling_quantile_rb : public rolling_rank_base_rb<T> {
    double percent{0};
    explicit rolling_quantile_rb(int size, double percent_) : rolling_rank_base_rb<T>(size), percent{percent_} {}

    using rolling_rank_base_rb<T>::m_valid_count;
    using rolling_rank_base_rb<T>::handle;
    using rolling_rank_base_rb<T>::m_sorted_data;

    double operator()(T new_value) {
        handle(new_value);
        if (m_valid_count > 0) {
            double idx = (m_valid_count - 1) * percent;
            long nth_lb = std::lround(std::floor(idx));
            long nth_ub = std::lround(std::ceil(idx));
            if (nth_lb < 0) nth_lb = 0;
            if (nth_lb >= m_valid_count) nth_lb = m_valid_count - 1;
            if (nth_ub < 0) nth_ub = 0;
            if (nth_ub >= m_valid_count) nth_ub = m_valid_count - 1;
            if (nth_lb == nth_ub) return m_sorted_data[nth_lb].data;
            return m_sorted_data[nth_lb].data * (nth_ub - idx) + m_sorted_data[nth_ub].data * (idx - nth_lb);
        } else
            return NAN;
    }
};

template <typename T>
struct rolling_rank_count_rb {
    int window_size;
    int m_count{0};
    std::vector<T> m_container;

    explicit rolling_rank_count_rb(int size) : window_size{size} { m_container.resize(window_size, T()); }

    double operator()(T new_value) {
        m_container[(m_count++) % window_size] = new_value;
        if (std::isnan(new_value) || m_count < window_size) {
            return NAN;
        }
        int nlte = 0, neq = 0, nv = 0;
        for (T val : m_container) {
            if (std::isnan(val)) continue;
            if (val <= new_value) nlte++;
            if (val == new_value) neq++;
            nv++;
        }
        // return nv <= 1 ? NAN : (nlte - 1.0) / (nv - 1.0);
        return nv <= 1 ? NAN : (2 * nlte - neq - 1.0) / (2.0 * (nv - 1.0));
    }
};

struct rolling_regression2_rb {
    double sum_x{0}, sum_y{0}, sum_x2{0}, sum_xy{0};
    double a{NAN}, b{NAN};
    int window_size;
    int m_count{0}, m_valid_count{0};
    int m_head_index{0};
    std::vector<std::pair<double, double>> m_container;

    rolling_regression2_rb(int size) : window_size{size + 1} {
        m_container.resize(window_size, std::pair<double, double>());
    }

    void delete_old() {
        int old_index = m_head_index - window_size;
        if (old_index < 0) old_index += window_size;
        const auto& old_value = m_container[old_index];
        if (std::isfinite(old_value.first) && std::isfinite(old_value.second)) {
            sum_x -= old_value.second;
            sum_y -= old_value.first;
            sum_x2 -= old_value.second * old_value.second;
            sum_xy -= old_value.first * old_value.second;
            --m_valid_count;
        }
    }
    void add_new() {
        const auto& new_value = m_container[m_head_index - 1];
        if (std::isfinite(new_value.first) && std::isfinite(new_value.second)) {
            sum_x += new_value.second;
            sum_y += new_value.first;
            sum_x2 += new_value.second * new_value.second;
            sum_xy += new_value.first * new_value.second;
            ++m_valid_count;
        }

        if (m_valid_count > 1) {
            b = (m_valid_count * sum_xy - sum_x * sum_y) / (m_valid_count * sum_x2 - sum_x * sum_x);
            a = (sum_y - b * sum_x) / m_valid_count;
        } else {
            a = NAN;
            b = NAN;
        }
    }
    void operator()(double y, double x) {
        m_container[m_head_index++] = {y, x};
        ++m_count;

        if (m_count >= window_size) {
            delete_old();
        }
        add_new();
        if (m_head_index == window_size) m_head_index = 0;
    }
};

struct rolling_regression3_rb {
    struct Data {
        double y{NAN}, x1{NAN}, x2{NAN};
        Data() = default;
        Data(double y_, double x1_, double x2_) : y{y_}, x1{x1_}, x2{x2_} {}
    };

    double sum_x1_2{0}, sum_x2_2{0}, sum_y_2{0};
    double sum_x1{0}, sum_x2{0}, sum_y{0};
    double sum_x12{0}, sum_x1y{0}, sum_x2y{0};  // cross term
    double b0{NAN}, b1{NAN}, b2{NAN};
    int window_size;
    int m_count{0}, m_valid_count{0};
    int m_head_index{0};
    std::vector<Data> m_container;

    rolling_regression3_rb(int size) : window_size{size + 1} { m_container.resize(window_size, Data()); }

    void delete_old() {
        int old_index = m_head_index - window_size;
        if (old_index < 0) old_index += window_size;
        const auto& old_value = m_container[old_index];
        if (std::isfinite(old_value.x1) && std::isfinite(old_value.x2) && std::isfinite(old_value.y)) {
            sum_x12 -= old_value.x1 * old_value.x2;
            sum_x1_2 -= old_value.x1 * old_value.x1;
            sum_x2_2 -= old_value.x2 * old_value.x2;
            sum_y_2 -= old_value.y * old_value.y;
            sum_x1y -= old_value.x1 * old_value.y;
            sum_x2y -= old_value.x2 * old_value.y;
            sum_x1 -= old_value.x1;
            sum_x2 -= old_value.x2;
            sum_y -= old_value.y;
            --m_valid_count;
        }
    }
    void add_new() {
        const auto& new_value = m_container[m_head_index - 1];
        if (std::isfinite(new_value.x1) && std::isfinite(new_value.x2) && std::isfinite(new_value.y)) {
            sum_x12 += new_value.x1 * new_value.x2;
            sum_x1_2 += new_value.x1 * new_value.x1;
            sum_x2_2 += new_value.x2 * new_value.x2;
            sum_y_2 += new_value.y * new_value.y;
            sum_x1y += new_value.x1 * new_value.y;
            sum_x2y += new_value.x2 * new_value.y;
            sum_x1 += new_value.x1;
            sum_x2 += new_value.x2;
            sum_y += new_value.y;
            ++m_valid_count;
        }

        if (m_valid_count < 3) {
            b0 = NAN;
            b1 = NAN;
            b2 = NAN;
        } else {
            double sum_X1_2 = sum_x1_2 - sum_x1 * sum_x1 / m_valid_count;
            double sum_X2_2 = sum_x2_2 - sum_x2 * sum_x2 / m_valid_count;
            double sum_X1Y = sum_x1y - sum_x1 * sum_y / m_valid_count;
            double sum_X2Y = sum_x2y - sum_x2 * sum_y / m_valid_count;
            double sum_X12 = sum_x12 - sum_x1 * sum_x2 / m_valid_count;

            double denominator = sum_X1_2 * sum_X2_2 - sum_X12 * sum_X12;
            b1 = (sum_X2_2 * sum_X1Y - sum_X12 * sum_X2Y) / denominator;
            b2 = (sum_X1_2 * sum_X2Y - sum_X12 * sum_X1Y) / denominator;
            b0 = (sum_y - b1 * sum_x1 - b2 * sum_x2) / m_valid_count;
        }
    }
    void operator()(double y, double x1, double x2) {
        m_container[m_head_index++] = {y, x1, x2};
        ++m_count;

        if (m_count >= window_size) {
            delete_old();
        }
        add_new();
        if (m_head_index == window_size) m_head_index = 0;
    }
};

struct rolling_ema_hl_rb : public rolling_rb_base<double> {
    double total_{0}, total_w{0}, result{0};
    double decay_coeff{0};

    rolling_ema_hl_rb(int size, double hl) : rolling_rb_base<double>(size) { decay_coeff = ema_hl2decay(hl); }

    void delete_old() {
        int old_index = get_old_index();
        const double& old_value = m_container[old_index];
        if (std::isfinite(old_value)) {
            total_ -= old_value * std::pow(decay_coeff, window_size - 2);
            total_w -= std::pow(decay_coeff, window_size - 2);
            --m_valid_count;
        }
    }
    void add_new() {
        const double& new_value = m_container[m_head_index - 1];
        if (std::isfinite(new_value)) {
            total_ = total_ * decay_coeff + new_value;
            total_w = total_w * decay_coeff + 1.0;
            ++m_valid_count;
        } else {
            total_ = total_ * decay_coeff;
            total_w = total_w * decay_coeff;
        }

        if (m_count > window_size - 2 && m_valid_count > 0)
            result = total_ / total_w;
        else
            result = NAN;
    }
    double operator()(double data) {
        m_container[m_head_index++] = data;
        ++m_count;

        if (m_count >= window_size) {
            delete_old();
        }
        add_new();
        if (m_head_index == window_size) m_head_index = 0;
        return result;
    }
};

}  // namespace ornate

#endif
