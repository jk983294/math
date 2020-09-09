#ifndef ORNATE_MATH_STATS_ROLLING_RB_RANGE_H
#define ORNATE_MATH_STATS_ROLLING_RB_RANGE_H

#include <algorithm>
#include <functional>
#include "math_container.h"

namespace ornate {

struct rolling_mean_rb_range {
    struct stat {
        long double total_sum{0};
        int m_valid_count{0};
        void clear() {
            total_sum = 0;
            m_valid_count = 0;
        }
        double calc() const {
            if (m_valid_count > 0)
                return total_sum / m_valid_count;
            else
                return NAN;
        }
    };
    int m_column_size;
    std::vector<stat> stats;

    explicit rolling_mean_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    template <typename T, typename TOut>
    void operator()(const T* old_row, const T* new_row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            auto& stat = stats[i];
            if (old_row) {
                auto old_value = old_row[i];
                if (std::isfinite(old_value)) {
                    stat.total_sum -= old_value;
                    --stat.m_valid_count;
                }
            }

            auto new_value = new_row[i];
            if (std::isfinite(new_value)) {
                stat.total_sum += new_value;
                ++stat.m_valid_count;
            }

            output[i] = stat.calc();
        }
    }

    void init() {
        for (auto& stat : stats) stat.clear();
    }

    template <typename T, typename TOut>
    void full_single(int idx, const T* _row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            if (std::isfinite(_row[i])) {
                ++stats[i].m_valid_count;
                stats[i].total_sum += _row[i];
            }
        }
    }

    template <typename TOut>
    void final_result(TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            output[i] = stats[i].calc();
        }
    }

    void set_row_size(int row) {}
    void set_param(const std::string& key, const std::string& value) {}
};

struct rolling_sum_rb_range {
    struct stat {
        long double total_sum{0};
        int m_valid_count{0};
        void clear() {
            total_sum = 0;
            m_valid_count = 0;
        }
        double calc() const {
            if (m_valid_count > 0)
                return total_sum;
            else
                return NAN;
        }
    };
    int m_column_size;
    std::vector<stat> stats;
    int sign{0};  // 0: all, 1: + only, -1: - only

    explicit rolling_sum_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    template <typename T, typename TOut>
    void operator()(const T* old_row, const T* new_row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            auto& stat = stats[i];
            if (old_row) {
                auto old_value = old_row[i];
                if (std::isfinite(old_value)) {
                    if (is_same_sign(old_value, sign)) {
                        stat.total_sum -= old_value;
                        --stat.m_valid_count;
                    }
                }
            }

            auto new_value = new_row[i];
            if (std::isfinite(new_value)) {
                if (is_same_sign(new_value, sign)) {
                    stat.total_sum += new_value;
                    ++stat.m_valid_count;
                }
            }

            if (stat.m_valid_count > 0)
                output[i] = stat.total_sum;
            else
                output[i] = NAN;
        }
    }

    void init() {
        for (auto& stat : stats) stat.clear();
    }

    template <typename T, typename TOut>
    void full_single(int idx, const T* _row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            if (std::isfinite(_row[i])) {
                if (is_same_sign(_row[i], sign)) {
                    ++stats[i].m_valid_count;
                    stats[i].total_sum += _row[i];
                }
            }
        }
    }

    template <typename TOut>
    void final_result(TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            output[i] = stats[i].calc();
        }
    }

    void set_row_size(int row) {}
    void set_param(const std::string& key, const std::string& value) {
        if (key == "sign") {
            sign = std::stoi(value);
        }
    }
};

struct rolling_prod_rb_range {
    struct stat {
        long double total_prod{1};
        int m_valid_count{0};
        void clear() {
            total_prod = 1;
            m_valid_count = 0;
        }
        double calc() const {
            if (m_valid_count > 0)
                return total_prod;
            else
                return NAN;
        }
    };
    int m_column_size;
    std::vector<stat> stats;

    explicit rolling_prod_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    template <typename T, typename TOut>
    void operator()(const T* old_row, const T* new_row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            auto& stat = stats[i];
            if (old_row) {
                auto old_value = old_row[i];
                if (std::isfinite(old_value)) {
                    stat.total_prod /= old_value;
                    --stat.m_valid_count;
                }
            }

            auto new_value = new_row[i];
            if (std::isfinite(new_value)) {
                stat.total_prod *= new_value;
                ++stat.m_valid_count;
            }

            output[i] = stat.calc();
        }
    }

    void init() {
        for (auto& stat : stats) stat.clear();
    }

    template <typename T, typename TOut>
    void full_single(int idx, const T* _row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            if (std::isfinite(_row[i])) {
                ++stats[i].m_valid_count;
                stats[i].total_prod *= _row[i];
            }
        }
    }

    template <typename TOut>
    void final_result(TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            output[i] = stats[i].calc();
        }
    }

    void set_row_size(int row) {}
    void set_param(const std::string& key, const std::string& value) {}
};

struct rolling_delay_rb_range {
    int m_column_size;

    explicit rolling_delay_rb_range(int column_size_) : m_column_size{column_size_} {}

    template <typename T, typename TOut>
    void operator()(const T* old_row, const T* new_row, TOut* output) {
        if (old_row) {
            detail::_data_copy2vector(old_row, output, m_column_size);
        } else {
            std::fill(output, output + m_column_size, NAN);
        }
    }

    void set_row_size(int row) {}
    void set_param(const std::string& key, const std::string& value) {}
};

struct rolling_delta_rb_range {
    int m_column_size;

    explicit rolling_delta_rb_range(int column_size_) : m_column_size{column_size_} {}

    template <typename T, typename TOut>
    void operator()(const T* old_row, const T* new_row, TOut* output) {
        if (old_row) {
            for (int ii = 0; ii < m_column_size; ++ii) {
                output[ii] = new_row[ii] - old_row[ii];
            }
        } else {
            std::fill(output, output + m_column_size, NAN);
        }
    }

    void set_row_size(int row) {}
    void set_param(const std::string& key, const std::string& value) {}
};

struct rolling_pct_rb_range {
    int m_column_size;

    explicit rolling_pct_rb_range(int column_size_) : m_column_size{column_size_} {}

    template <typename T, typename TOut>
    void operator()(const T* old_row, const T* new_row, TOut* output) {
        if (old_row) {
            for (int ii = 0; ii < m_column_size; ++ii) {
                output[ii] = new_row[ii] / old_row[ii] - 1.0;
            }
        } else {
            std::fill(output, output + m_column_size, NAN);
        }
    }

    void set_row_size(int row) {}
    void set_param(const std::string& key, const std::string& value) {}
};

struct rolling_variance_rb_range {
    struct stat {
        long double total_sum{0}, total_square_sum{0};
        int m_valid_count{0};

        void clear() {
            total_sum = 0;
            total_square_sum = 0;
            m_valid_count = 0;
        }
        double calc(bool demean_) const {
            if (m_valid_count > 1) {
                if (demean_) {
                    double mean = total_sum / m_valid_count;
                    return (total_square_sum - mean * mean * m_valid_count) / (m_valid_count - 1);
                } else {
                    return (total_square_sum) / m_valid_count;
                }
            }
            return NAN;
        }
    };
    int m_column_size;
    std::vector<stat> stats;
    bool demean{true};  // for some series, its mean is 0, like return, no need to - mean

    explicit rolling_variance_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    template <typename T, typename TOut>
    void operator()(const T* old_row, const T* new_row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            auto& st = stats[i];
            if (old_row) {
                delete_old(st, old_row[i]);
            }
            output[i] = add_new(st, new_row[i]);
        }
    }

    template <typename T>
    void delete_old(stat& st, T old_value) {
        if (std::isfinite(old_value)) {
            st.total_sum -= old_value;
            st.total_square_sum -= old_value * old_value;
            --st.m_valid_count;
        }
    }

    template <typename T>
    double add_new(stat& st, T new_value) {
        if (std::isfinite(new_value)) {
            st.total_sum += new_value;
            st.total_square_sum += new_value * new_value;
            ++st.m_valid_count;
        }
        return st.calc(demean);
    }

    void init() {
        for (auto& stat : stats) stat.clear();
    }

    template <typename T, typename TOut>
    void full_single(int idx, const T* _row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            if (std::isfinite(_row[i])) {
                ++stats[i].m_valid_count;
                stats[i].total_sum += _row[i];
                stats[i].total_square_sum += _row[i] * _row[i];
            }
        }
    }

    template <typename TOut>
    void final_result(TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            output[i] = stats[i].calc(demean);
        }
    }

    void set_row_size(int row) {}
    void set_param(const std::string& key, const std::string& value) {
        if (key == "demean") {
            demean = value == "true";
        }
    }
};

struct rolling_std_rb_range {
    struct stat {
        long double total_sum{0}, total_square_sum{0};
        int m_valid_count{0};

        void clear() {
            total_sum = total_square_sum = 0;
            m_valid_count = 0;
        }
        double calc(bool demean_) const {
            if (m_valid_count > 1) {
                if (demean_) {
                    double var = (total_square_sum - total_sum * total_sum / (m_valid_count)) / (m_valid_count - 1);
                    return std::sqrt(var);
                } else {
                    return std::sqrt((total_square_sum) / m_valid_count);
                }
            }
            return NAN;
        }
    };
    int m_column_size;
    std::vector<stat> stats;
    bool demean{true};  // for some series, its mean is 0, like return, no need to - mean

    explicit rolling_std_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    template <typename T, typename TOut>
    void operator()(const T* old_row, const T* new_row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            auto& st = stats[i];
            if (old_row) {
                delete_old(st, old_row[i]);
            }
            output[i] = add_new(st, new_row[i]);
        }
    }

    template <typename T>
    void delete_old(stat& st, T old_value) {
        if (std::isfinite(old_value)) {
            st.total_sum -= old_value;
            st.total_square_sum -= old_value * old_value;
            --st.m_valid_count;
        }
    }

    template <typename T>
    double add_new(stat& st, T new_value) {
        if (std::isfinite(new_value)) {
            st.total_sum += new_value;
            st.total_square_sum += new_value * new_value;
            ++st.m_valid_count;
        }
        return st.calc(demean);
    }

    void init() {
        for (auto& stat : stats) stat.clear();
    }

    template <typename T, typename TOut>
    void full_single(int idx, const T* _row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            if (std::isfinite(_row[i])) {
                ++stats[i].m_valid_count;
                stats[i].total_sum += _row[i];
                stats[i].total_square_sum += _row[i] * _row[i];
            }
        }
    }

    template <typename TOut>
    void final_result(TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            output[i] = stats[i].calc(demean);
        }
    }

    void set_row_size(int row) {}
    void set_param(const std::string& key, const std::string& value) {
        if (key == "demean") {
            demean = value == "true";
        }
    }
};

struct rolling_zscore_rb_range {
    struct stat {
        long double total_sum{0}, total_square_sum{0};
        int m_valid_count{0};
        void clear() {
            total_sum = 0;
            total_square_sum = 0;
            m_valid_count = 0;
        }
        double calc(double latest_val) const {
            if (std::isfinite(latest_val)) {
                double mean = total_sum / m_valid_count;
                double stddev = std::sqrt((total_square_sum - mean * mean * m_valid_count) / (m_valid_count - 1));
                return (latest_val - mean) / stddev;
            }
            return NAN;
        }
    };
    int m_column_size;
    std::vector<stat> stats;

    explicit rolling_zscore_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    template <typename T, typename TOut>
    void operator()(const T* old_row, const T* new_row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            auto& st = stats[i];
            if (old_row) {
                delete_old(st, old_row[i]);
            }
            output[i] = add_new(st, new_row[i]);
        }
    }

    template <typename T>
    void delete_old(stat& st, T old_value) {
        if (std::isfinite(old_value)) {
            st.total_sum -= old_value;
            st.total_square_sum -= old_value * old_value;
            --st.m_valid_count;
        }
    }

    template <typename T>
    double add_new(stat& st, T new_value) {
        if (std::isfinite(new_value)) {
            st.total_sum += new_value;
            st.total_square_sum += new_value * new_value;
            ++st.m_valid_count;
            double mean = st.total_sum / st.m_valid_count;
            double stddev = std::sqrt((st.total_square_sum - mean * mean * st.m_valid_count) / (st.m_valid_count - 1));
            return (new_value - mean) / stddev;
        } else {
            return NAN;
        }
    }

    void init() {
        for (auto& stat : stats) stat.clear();
    }

    template <typename T, typename TOut>
    void full_single(int idx, const T* _row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            if (std::isfinite(_row[i])) {
                ++stats[i].m_valid_count;
                stats[i].total_sum += _row[i];
                stats[i].total_square_sum += _row[i] * _row[i];
            }
        }
        if (idx == 0) {
            detail::_data_copy2vector(_row, output, m_column_size);
        }
    }

    template <typename TOut>
    void final_result(TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            output[i] = stats[i].calc(output[i]);
        }
    }

    void set_row_size(int row) {}
    void set_param(const std::string& key, const std::string& value) {}
};

struct rolling_score_rb_range {
    struct stat {
        long double total_sum{0};
        int m_valid_count{0};
        void clear() {
            total_sum = 0;
            m_valid_count = 0;
        }
        double calc(double latest_val) const {
            if (std::isfinite(latest_val)) {
                return latest_val - total_sum / m_valid_count;
            }
            return NAN;
        }
    };
    int m_column_size;
    std::vector<stat> stats;

    explicit rolling_score_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    template <typename T, typename TOut>
    void operator()(const T* old_row, const T* new_row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            auto& stat = stats[i];
            if (old_row) {
                auto old_value = old_row[i];
                if (std::isfinite(old_value)) {
                    stat.total_sum -= old_value;
                    --stat.m_valid_count;
                }
            }

            auto new_value = new_row[i];
            if (std::isfinite(new_value)) {
                stat.total_sum += new_value;
                ++stat.m_valid_count;
                output[i] = new_value - stat.total_sum / stat.m_valid_count;
            } else
                output[i] = NAN;
        }
    }

    void init() {
        for (auto& stat : stats) stat.clear();
    }

    template <typename T, typename TOut>
    void full_single(int idx, const T* _row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            if (std::isfinite(_row[i])) {
                ++stats[i].m_valid_count;
                stats[i].total_sum += _row[i];
            }
        }
        if (idx == 0) {
            detail::_data_copy2vector(_row, output, m_column_size);
        }
    }

    template <typename TOut>
    void final_result(TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            output[i] = stats[i].calc(output[i]);
        }
    }

    void set_row_size(int row) {}
    void set_param(const std::string& key, const std::string& value) {}
};

struct rolling_cov_rb_range {
    struct stat {
        long double sumx{0}, sumy{0}, sumxy{0};
        int m_valid_count{0};
        void clear() {
            sumx = sumy = sumxy = 0;
            m_valid_count = 0;
        }
        double calc() const {
            if (m_valid_count > 1) {
                double mean_x = sumx / m_valid_count;
                double mean_y = sumy / m_valid_count;
                return (sumxy - mean_x * mean_y * m_valid_count) / (m_valid_count - 1);
            }
            return NAN;
        }
    };
    int m_column_size;
    std::vector<stat> stats;

    explicit rolling_cov_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    template <typename T, typename TOut>
    void operator()(const T* old_row0, const T* old_row1, const T* new_row0, const T* new_row1, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            auto& st = stats[i];
            if (old_row0) {
                delete_old(st, old_row0[i], old_row1[i]);
            }
            output[i] = add_new(st, new_row0[i], new_row1[i]);
        }
    }

    template <typename T>
    void delete_old(stat& st, T old_value0, T old_value1) {
        if (std::isfinite(old_value0) && std::isfinite(old_value1)) {
            st.sumxy -= old_value0 * old_value1;
            st.sumx -= old_value0;
            st.sumy -= old_value1;
            --st.m_valid_count;
        }
    }
    template <typename T>
    double add_new(stat& st, T data0, T data1) {
        if (std::isfinite(data0) && std::isfinite(data1)) {
            st.sumxy += data0 * data1;
            st.sumx += data0;
            st.sumy += data1;
            ++st.m_valid_count;
        }

        return st.calc();
    }

    void init() {
        for (auto& stat : stats) stat.clear();
    }

    template <typename T, typename TOut>
    void full_single(int idx, const T* x_row, const T* y_row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            if (std::isfinite(x_row[i]) && std::isfinite(y_row[i])) {
                ++stats[i].m_valid_count;
                stats[i].sumxy += x_row[i] * y_row[i];
                stats[i].sumx += x_row[i];
                stats[i].sumy += y_row[i];
            }
        }
    }

    template <typename TOut>
    void final_result(TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            output[i] = stats[i].calc();
        }
    }

    void set_row_size(int row) {}
    void set_param(const std::string& key, const std::string& value) {}
};

struct rolling_corr_rb_range {
    struct stat {
        long double sumx{0}, sumy{0}, sumxy{0}, sum_x2{0}, sum_y2{0};
        int m_valid_count{0};
        void clear() {
            sumx = sumy = sumxy = sum_x2 = sum_y2 = 0;
            m_valid_count = 0;
        }
        double calc() const {
            if (m_valid_count > 0) {
                double cov = (sumxy / m_valid_count - sumx * sumy / (m_valid_count * m_valid_count));
                double var1 = (sum_x2 / m_valid_count - sumx * sumx / (m_valid_count * m_valid_count));
                double var2 = (sum_y2 / m_valid_count - sumy * sumy / (m_valid_count * m_valid_count));
                double numerator = sqrt(var1 * var2);
                if (numerator > 0) return cov / numerator;
            }
            return NAN;
        }
    };
    int m_column_size;
    std::vector<stat> stats;

    explicit rolling_corr_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    template <typename T, typename TOut>
    void operator()(const T* old_row0, const T* old_row1, const T* new_row0, const T* new_row1, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            auto& st = stats[i];
            if (old_row0) {
                delete_old(st, old_row0[i], old_row1[i]);
            }
            output[i] = add_new(st, new_row0[i], new_row1[i]);
        }
    }

    template <typename T>
    void delete_old(stat& st, T old_value0, T old_value1) {
        if (std::isfinite(old_value0) && std::isfinite(old_value1)) {
            st.sumxy -= old_value0 * old_value1;
            st.sumx -= old_value0;
            st.sumy -= old_value1;
            st.sum_x2 -= old_value0 * old_value0;
            st.sum_y2 -= old_value1 * old_value1;
            --st.m_valid_count;
        }
    }

    template <typename T>
    double add_new(stat& st, T data0, T data1) {
        if (std::isfinite(data0) && std::isfinite(data1)) {
            st.sumxy += data0 * data1;
            st.sumx += data0;
            st.sumy += data1;
            st.sum_x2 += data0 * data0;
            st.sum_y2 += data1 * data1;
            ++st.m_valid_count;
        }
        return st.calc();
    }

    void init() {
        for (auto& stat : stats) stat.clear();
    }

    template <typename T, typename TOut>
    void full_single(int idx, const T* x_row, const T* y_row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            if (std::isfinite(x_row[i]) && std::isfinite(y_row[i])) {
                ++stats[i].m_valid_count;
                stats[i].sumxy += x_row[i] * y_row[i];
                stats[i].sumx += x_row[i];
                stats[i].sumy += y_row[i];
                stats[i].sum_x2 += x_row[i] * x_row[i];
                stats[i].sum_y2 += y_row[i] * y_row[i];
            }
        }
    }

    template <typename TOut>
    void final_result(TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            output[i] = stats[i].calc();
        }
    }

    void set_row_size(int row) {}
    void set_param(const std::string& key, const std::string& value) {}
};

struct rolling_skew_rb_range {
    struct stat {
        long double total_x1{0}, total_x2{0}, total_x3{0};
        int m_valid_count{0};
        void clear() {
            total_x1 = total_x2 = total_x3 = 0;
            m_valid_count = 0;
        }
        double calc() const {
            if (m_valid_count >= 2) {
                // long double mx = total_x1 / m_valid_count;
                // long double d = pow(total_x2 - m_valid_count * mx * mx, 1.5);
                // if(d > 0)
                //   return sqrt(m_valid_count) * (total_x3 - 3 * mx * total_x2 + 2 * m_valid_count * mx * mx * mx) / d;
                // else return NAN;
                long double mean = total_x1 / m_valid_count;
                long double mean2 = mean * mean;
                long double var = total_x2 / m_valid_count - mean2;
                if (var > 1e-14) {
                    long double mean3 = mean * mean * mean;
                    long double m3 = total_x3 / m_valid_count - 3 * mean * total_x2 / m_valid_count + 2 * mean3;
                    return m3 / std::pow(var, 1.5);
                }
            }
            return NAN;
        }
    };
    int m_column_size;
    std::vector<stat> stats;

    explicit rolling_skew_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    template <typename T, typename TOut>
    void operator()(const T* old_row, const T* new_row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            auto& st = stats[i];
            if (old_row) {
                delete_old(st, old_row[i]);
            }
            output[i] = add_new(st, new_row[i]);
        }
    }

    template <typename T>
    void delete_old(stat& st, T old_value) {
        if (std::isfinite(old_value)) {
            st.total_x1 -= old_value;
            st.total_x2 -= old_value * old_value;
            st.total_x3 -= old_value * old_value * old_value;
            --st.m_valid_count;
        }
    }

    template <typename T>
    double add_new(stat& st, T new_value) {
        if (std::isfinite(new_value)) {
            st.total_x1 += new_value;
            st.total_x2 += new_value * new_value;
            st.total_x3 += new_value * new_value * new_value;
            ++st.m_valid_count;
        }

        return st.calc();
    }

    void init() {
        for (auto& stat : stats) stat.clear();
    }

    template <typename T, typename TOut>
    void full_single(int idx, const T* _row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            if (std::isfinite(_row[i])) {
                ++stats[i].m_valid_count;
                stats[i].total_x1 += _row[i];
                stats[i].total_x2 += _row[i] * _row[i];
                stats[i].total_x3 += _row[i] * _row[i] * _row[i];
            }
        }
    }

    template <typename TOut>
    void final_result(TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            output[i] = stats[i].calc();
        }
    }

    void set_row_size(int row) {}
    void set_param(const std::string& key, const std::string& value) {}
};

struct rolling_kurtosis_rb_range {
    struct stat {
        long double total_x1{0}, total_x2{0}, total_x3{0}, total_x4{0};
        int m_valid_count{0};
        void clear() {
            total_x1 = total_x2 = total_x3 = total_x4 = 0;
            m_valid_count = 0;
        }
        double calc() const {
            if (m_valid_count >= 2) {
                double mean = total_x1 / m_valid_count;
                double mean2 = mean * mean;
                double var = total_x2 / m_valid_count - mean2;
                if (var > 1e-14) {
                    double mean3 = mean2 * mean;
                    double mean4 = mean3 * mean;
                    double m4 = total_x4 / m_valid_count - 4 * mean * total_x3 / m_valid_count +
                                6 * mean2 * total_x2 / m_valid_count - 3 * mean4;
                    return m4 / std::pow(var, 2) - 3.0;
                }
            }
            return NAN;
        }
    };
    int m_column_size;
    std::vector<stat> stats;

    explicit rolling_kurtosis_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    template <typename T, typename TOut>
    void operator()(const T* old_row, const T* new_row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            auto& st = stats[i];
            if (old_row) {
                delete_old(st, old_row[i]);
            }
            output[i] = add_new(st, new_row[i]);
        }
    }

    template <typename T>
    void delete_old(stat& st, T old_value) {
        if (std::isfinite(old_value)) {
            st.total_x1 -= old_value;
            st.total_x2 -= old_value * old_value;
            st.total_x3 -= old_value * old_value * old_value;
            st.total_x4 -= old_value * old_value * old_value * old_value;
            --st.m_valid_count;
        }
    }

    template <typename T>
    double add_new(stat& st, T new_value) {
        if (std::isfinite(new_value)) {
            st.total_x1 += new_value;
            st.total_x2 += new_value * new_value;
            st.total_x3 += new_value * new_value * new_value;
            st.total_x4 += new_value * new_value * new_value * new_value;
            ++st.m_valid_count;
        }

        return st.calc();
    }

    void init() {
        for (auto& stat : stats) stat.clear();
    }

    template <typename T, typename TOut>
    void full_single(int idx, const T* _row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            if (std::isfinite(_row[i])) {
                ++stats[i].m_valid_count;
                stats[i].total_x1 += _row[i];
                stats[i].total_x2 += _row[i] * _row[i];
                stats[i].total_x3 += _row[i] * _row[i] * _row[i];
                stats[i].total_x4 += _row[i] * _row[i] * _row[i] * _row[i];
            }
        }
    }

    template <typename TOut>
    void final_result(TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            output[i] = stats[i].calc();
        }
    }

    void set_row_size(int row) {}
    void set_param(const std::string& key, const std::string& value) {}
};

struct rolling_decay_rb_range {
    struct stat {
        long double total_x1{0}, total{0};
        int m_valid_count{0}, m_valid_x1_count{0};
        void clear() {
            total = total_x1 = 0;
            m_valid_count = m_valid_x1_count = 0;
        }
        double calc() const {
            if (m_valid_count > 0)
                return total / m_valid_count;
            else
                return NAN;
        }
    };
    int m_column_size;
    int m_row_size{0};  // window size
    int count{0};
    std::vector<stat> stats;

    explicit rolling_decay_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    void set_row_size(int row) { m_row_size = row; }
    void set_param(const std::string& key, const std::string& value) {}

    template <typename T, typename TOut>
    void operator()(const T* old_row, const T* new_row, TOut* output) {
        ++count;
        for (int i = 0; i < m_column_size; ++i) {
            auto& st = stats[i];
            if (old_row) {
                delete_old(st, old_row[i]);
            }
            output[i] = add_new(st, new_row[i]);
        }
    }

    template <typename T>
    void delete_old(stat& st, T old_value) {
        if (std::isfinite(old_value)) {
            st.total -= old_value;
            st.total_x1 -= old_value;
            --st.m_valid_count;
            --st.m_valid_x1_count;
        }
    }

    template <typename T>
    double add_new(stat& st, T new_value) {
        if (count > m_row_size) {
            st.total -= st.total_x1;
            st.m_valid_count -= st.m_valid_x1_count;
        }
        if (std::isfinite(new_value)) {
            st.total_x1 += new_value;
            ++st.m_valid_x1_count;
            if (count <= m_row_size) {
                st.total += new_value * count;
                st.m_valid_count += count;
            } else {
                st.total += new_value * m_row_size;
                st.m_valid_count += m_row_size;
            }
        }
        return st.calc();
    }

    void init() {
        count = 0;
        for (auto& stat : stats) stat.clear();
    }

    template <typename T, typename TOut>
    void full_single(int idx, const T* _row, TOut* output) {
        ++count;
        for (int i = 0; i < m_column_size; ++i) {
            if (std::isfinite(_row[i])) {
                stats[i].m_valid_count += (m_row_size - idx);
                stats[i].total += _row[i] * (m_row_size - idx);
                stats[i].total_x1 += _row[i];
                ++stats[i].m_valid_x1_count;
            }
        }
    }

    template <typename TOut>
    void final_result(TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            output[i] = stats[i].calc();
        }
    }
};

template <typename TData, typename TCmp = std::greater<TData>>
struct rolling_mq_rb_range {
public:
    struct Cell {
        int seq{-1};
        TData data;

        Cell() : data{TData()} {}
    };

    struct stat {
        int front = 0;
        int rear = 0;
        int seq = -1;
    };
    int capacity{0};
    int m_column_size;
    std::vector<stat> stats;
    std::vector<Cell> _data;
    TCmp cmp;

    explicit rolling_mq_rb_range(int column_size_) : m_column_size{column_size_} {
        stats.resize(m_column_size);
        cmp = TCmp();
    }
    void set_row_size(int row) {
        capacity = row + 1;
        _data.resize(capacity * m_column_size);
    }
    void set_param(const std::string& key, const std::string& value) {}

    void push(TData value, stat& st, Cell* start_cell) {
        int oldest_seq = st.seq - capacity + 2;
        int ptr = st.front;
        for (; ptr != st.rear && start_cell[ptr].seq <= oldest_seq; ptr = (ptr + 1) % capacity)
            ;
        st.front = ptr;

        if (!isvalid(value)) {
            ++st.seq;
            return;
        }

        ptr = (st.rear + capacity - 1) % capacity;
        int end_ptr = (st.front + capacity - 1) % capacity;
        for (; ptr != end_ptr && cmp(start_cell[ptr].data, value); ptr = (ptr + capacity - 1) % capacity)
            ;
        ptr = (ptr + 1) % capacity;

        auto& cell = start_cell[ptr];
        cell.seq = ++st.seq;
        cell.data = value;
        st.rear = (ptr + 1) % capacity;
    }

    double top_index(stat& st, Cell* start_cell) {
        if (st.front == st.rear) return NAN;
        return st.seq - start_cell[st.front].seq;
    }

    double top_percent(stat& st, Cell* start_cell) {
        if (st.front == st.rear) return NAN;
        return double(st.seq - start_cell[st.front].seq) / (capacity - 2);
    }

    TData top(stat& st, Cell* start_cell) {
        if (st.front == st.rear) return get_nan<TData>();
        return start_cell[st.front].data;
    }
};

template <typename TData, typename TCmp = std::greater<TData>>
struct rolling_mq_index_rb_range : public rolling_mq_rb_range<TData, TCmp> {
    using rolling_mq_rb_range<TData, TCmp>::m_column_size;
    using rolling_mq_rb_range<TData, TCmp>::stats;
    using rolling_mq_rb_range<TData, TCmp>::Cell;
    using rolling_mq_rb_range<TData, TCmp>::_data;
    using rolling_mq_rb_range<TData, TCmp>::capacity;
    using rolling_mq_rb_range<TData, TCmp>::push;
    using rolling_mq_rb_range<TData, TCmp>::top_index;

    explicit rolling_mq_index_rb_range(int column_size_) : rolling_mq_rb_range<TData, TCmp>(column_size_) {}

    template <typename T, typename TOut>
    void operator()(const T* old_row, const T* new_row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            auto& st = stats[i];
            auto* start_cell = _data.data() + i * capacity;
            push(new_row[i], st, start_cell);
            output[i] = top_index(st, start_cell);
        }
    }
};

template <typename TData, typename TCmp = std::greater<TData>>
struct rolling_mq_percent_rb_range : public rolling_mq_rb_range<TData, TCmp> {
    using rolling_mq_rb_range<TData, TCmp>::m_column_size;
    using rolling_mq_rb_range<TData, TCmp>::stats;
    using rolling_mq_rb_range<TData, TCmp>::Cell;
    using rolling_mq_rb_range<TData, TCmp>::_data;
    using rolling_mq_rb_range<TData, TCmp>::capacity;
    using rolling_mq_rb_range<TData, TCmp>::push;
    using rolling_mq_rb_range<TData, TCmp>::top_percent;

    explicit rolling_mq_percent_rb_range(int column_size_) : rolling_mq_rb_range<TData, TCmp>(column_size_) {}

    template <typename T, typename TOut>
    void operator()(const T* old_row, const T* new_row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            auto& st = stats[i];
            auto* start_cell = _data.data() + i * capacity;
            push(new_row[i], st, start_cell);
            output[i] = top_percent(st, start_cell);
        }
    }
};

template <typename TData, typename TCmp = std::greater<TData>>
struct rolling_mq_value_rb_range : public rolling_mq_rb_range<TData, TCmp> {
    using rolling_mq_rb_range<TData, TCmp>::m_column_size;
    using rolling_mq_rb_range<TData, TCmp>::stats;
    using rolling_mq_rb_range<TData, TCmp>::Cell;
    using rolling_mq_rb_range<TData, TCmp>::_data;
    using rolling_mq_rb_range<TData, TCmp>::capacity;
    using rolling_mq_rb_range<TData, TCmp>::push;
    using rolling_mq_rb_range<TData, TCmp>::top;

    explicit rolling_mq_value_rb_range(int column_size_) : rolling_mq_rb_range<TData, TCmp>(column_size_) {}

    template <typename T, typename TOut>
    void operator()(const T* old_row, const T* new_row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            auto& st = stats[i];
            auto* start_cell = _data.data() + i * capacity;
            push(new_row[i], st, start_cell);
            output[i] = top(st, start_cell);
        }
    }
};

template <typename T>
struct rolling_rank_base_rb_range {
public:
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

    struct stat {
        int m_valid_count = 0;
    };
    int window_size{0}, m_count{0};
    int m_column_size;
    std::vector<stat> stats;
    std::vector<SortItem> m_sorted_data_;

    explicit rolling_rank_base_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }
    void set_row_size(int row) {
        window_size = row + 1;
        m_sorted_data_.resize(row * m_column_size);
    }

    std::tuple<int, int, int> handle(T new_value, T old_value, stat& st, SortItem* start_sorted_data) {
        int idx = -1, lower = -1;
        if (m_count >= window_size) {
            if (isvalid(old_value)) {
                int insert_pos_index = find_insert_point(start_sorted_data, st, old_value, m_count - window_size);

                if (isvalid(new_value)) {
                    if (new_value > old_value) {
                        idx = insert_right(start_sorted_data, st, new_value, insert_pos_index);
                    } else {
                        idx = insert_left(start_sorted_data, st, new_value, insert_pos_index);
                    }

                    lower = lower_bound(start_sorted_data, new_value, idx);
                } else {
                    std::rotate(start_sorted_data + insert_pos_index, start_sorted_data + insert_pos_index + 1,
                                start_sorted_data + st.m_valid_count);
                    --st.m_valid_count;
                }
            } else {
                if (isvalid(new_value)) {
                    idx = insert_left(start_sorted_data, st, new_value, st.m_valid_count);
                    lower = lower_bound(start_sorted_data, new_value, idx);
                    st.m_valid_count++;
                }
            }
        } else {
            if (isvalid(new_value)) {
                idx = insert_left(start_sorted_data, st, new_value, st.m_valid_count);
                lower = lower_bound(start_sorted_data, new_value, idx);
                st.m_valid_count++;
            }
        }

        return {lower, idx, st.m_valid_count};
    }

    int lower_bound(SortItem* start_sorted_data, T data, int idx) {
        while (idx >= 0 && start_sorted_data[idx].data == data) --idx;
        return idx + 1;
    }

    int insert_left(SortItem* start_sorted_data, stat& st, T data, int idx) {
        start_sorted_data[idx].data = data;
        start_sorted_data[idx].seq = m_count - 1;
        auto itr = std::upper_bound(start_sorted_data, start_sorted_data + idx, start_sorted_data[idx]);
        std::rotate(itr, start_sorted_data + idx, start_sorted_data + idx + 1);
        return itr - start_sorted_data;
    }

    int insert_right(SortItem* start_sorted_data, stat& st, T data, int idx) {
        start_sorted_data[idx].data = data;
        start_sorted_data[idx].seq = m_count - 1;
        auto itr =
            std::upper_bound(start_sorted_data + idx, start_sorted_data + st.m_valid_count, start_sorted_data[idx]);
        std::rotate(start_sorted_data + idx, start_sorted_data + idx + 1, itr);
        return int(itr - start_sorted_data) - 1;
    }

    int find_insert_point(SortItem* start_sorted_data, stat& st, T data, int seq) {
        auto itr = std::lower_bound(start_sorted_data, start_sorted_data + st.m_valid_count, data);
        int pos = itr - start_sorted_data;
        while (pos < st.m_valid_count && start_sorted_data[pos].seq != seq && start_sorted_data[pos].data == data)
            pos++;
        return pos;
    }
};

template <typename T>
struct rolling_quantile_rb_range : public rolling_rank_base_rb_range<T> {
    using rolling_rank_base_rb_range<T>::m_count;
    using rolling_rank_base_rb_range<T>::m_column_size;
    using rolling_rank_base_rb_range<T>::stats;
    using rolling_rank_base_rb_range<T>::m_sorted_data_;
    using rolling_rank_base_rb_range<T>::window_size;
    using rolling_rank_base_rb_range<T>::handle;
    double percent{0};

    explicit rolling_quantile_rb_range(int size, double percent_)
        : rolling_rank_base_rb_range<T>(size), percent{percent_} {}
    explicit rolling_quantile_rb_range(int size) : rolling_rank_base_rb_range<T>(size) {}

    template <typename TOut>
    void operator()(const T* old_row, const T* new_row, TOut* output) {
        ++m_count;
        for (int i = 0; i < m_column_size; ++i) {
            auto& st = stats[i];
            auto* start_sorted_data = m_sorted_data_.data() + i * (window_size - 1);

            T old_value = get_nan<T>();
            if (old_row) old_value = old_row[i];
            int _valid_count;
            std::tie(std::ignore, std::ignore, _valid_count) = handle(new_row[i], old_value, st, start_sorted_data);

            if (_valid_count > 0) {
                double idx = (_valid_count - 1) * percent;
                long nth_lb = std::lround(std::floor(idx));
                long nth_ub = std::lround(std::ceil(idx));
                if (nth_lb == nth_ub)
                    output[i] = start_sorted_data[nth_lb].data;
                else
                    output[i] = start_sorted_data[nth_lb].data * (nth_ub - idx) +
                                start_sorted_data[nth_ub].data * (idx - nth_lb);
            } else
                output[i] = NAN;
        }
    }

    void set_param(const std::string& key, const std::string& value) {
        if (key == "percent" || key == "arg1") {
            percent = std::stod(value);
        }
    }
};

template <typename T>
struct rolling_rank_count_rb_range {
    int window_size{0}, m_count{0};
    int m_column_size{0};
    std::vector<T> m_container;

    explicit rolling_rank_count_rb_range(int column_size_) : m_column_size{column_size_} {}
    void set_row_size(int row) {
        window_size = row;
        m_container.resize(window_size * m_column_size);
    }

    double calc(const T* x, T new_value) {
        if (std::isnan(new_value) || m_count < window_size) {
            return NAN;
        }
        int nlte = 0, neq = 0, nv = 0;
        for (int i = 0; i < window_size; ++i) {
            const T val = x[i];
            if (std::isnan(val)) continue;
            if (val <= new_value) nlte++;
            if (val == new_value) neq++;
            nv++;
        }
        // return nv <= 1 ? NAN : (nlte - 1.0) / (nv - 1.0);
        return nv <= 1 ? NAN : (2 * nlte - neq - 1.0) / (2.0 * (nv - 1.0));
    }

    template <typename TOut>
    void operator()(const T* old_row, const T* new_row, TOut* output) {
        int pos = m_count % window_size;
        ++m_count;
        for (int i = 0; i < m_column_size; ++i) {
            auto* start_data = m_container.data() + i * window_size;
            const T new_val = new_row[i];
            start_data[pos] = new_val;
            output[i] = calc(start_data, new_val);
        }
    }

    void set_param(const std::string& key, const std::string& value) {}
};

struct rolling_regression2_rb_range {
    struct stat {
        long double sum_x{0}, sum_y{0}, sum_x2{0}, sum_xy{0};
        int m_valid_count{0};
        void clear() {
            sum_x = sum_y = sum_x2 = sum_xy = 0;
            m_valid_count = 0;
        }
        std::pair<double, double> calc() const {
            double a = NAN, b = NAN;
            if (m_valid_count > 1) {
                b = (m_valid_count * sum_xy - sum_x * sum_y) / (m_valid_count * sum_x2 - sum_x * sum_x);
                a = (sum_y - b * sum_x) / m_valid_count;
            }
            return {a, b};
        }
    };
    int m_column_size;
    std::vector<stat> stats;

    explicit rolling_regression2_rb_range(int column_size_) : m_column_size{column_size_} {
        stats.resize(m_column_size);
    }

    template <typename T, typename TOut>
    void operator()(const T* old_y, const T* old_x, const T* new_y, const T* new_x, TOut* a, TOut* b) {
        for (int i = 0; i < m_column_size; ++i) {
            auto& st = stats[i];
            if (old_y) {
                delete_old(st, old_y[i], old_x[i]);
            }
            add_new(st, new_y[i], new_x[i]);
            std::tie(a[i], b[i]) = st.calc();
        }
    }

    template <typename T>
    void delete_old(stat& st, T old_y, T old_x) {
        if (std::isfinite(old_y) && std::isfinite(old_x)) {
            st.sum_x -= old_x;
            st.sum_y -= old_y;
            st.sum_x2 -= old_x * old_x;
            st.sum_xy -= old_y * old_x;
            --st.m_valid_count;
        }
    }
    template <typename T>
    void add_new(stat& st, T y, T x) {
        if (std::isfinite(y) && std::isfinite(x)) {
            st.sum_x += x;
            st.sum_y += y;
            st.sum_x2 += x * x;
            st.sum_xy += x * y;
            ++st.m_valid_count;
        }
    }

    void init() {
        for (auto& stat : stats) stat.clear();
    }

    template <typename T>
    void full_single(int idx, const T* y_row, const T* x_row) {
        for (int i = 0; i < m_column_size; ++i) {
            auto& st = stats[i];
            add_new(st, y_row[i], x_row[i]);
        }
    }

    template <typename TOut>
    void final_result(TOut* a, TOut* b) {
        for (int i = 0; i < m_column_size; ++i) {
            auto& st = stats[i];
            std::tie(a[i], b[i]) = st.calc();
        }
    }

    void set_row_size(int row) {}
    void set_param(const std::string& key, const std::string& value) {}
};

struct rolling_regression3_rb_range {
    struct stat {
        long double sum_x1_2{0}, sum_x2_2{0}, sum_y_2{0};
        long double sum_x1{0}, sum_x2{0}, sum_y{0};
        long double sum_x12{0}, sum_x1y{0}, sum_x2y{0};  // cross term
        // double b0{NAN}, b1{NAN}, b2{NAN};
        int m_valid_count{0};

        void clear() {
            sum_x1_2 = sum_x2_2 = sum_y_2 = sum_x1 = sum_x2 = sum_y = sum_x12 = sum_x1y = sum_x2y = 0;
            m_valid_count = 0;
        }
        std::tuple<double, double, double> calc() const {
            double b0 = NAN, b1 = NAN, b2 = NAN;
            if (m_valid_count >= 3) {
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
            return {b0, b1, b2};
        }
    };
    int m_column_size;
    std::vector<stat> stats;

    explicit rolling_regression3_rb_range(int column_size_) : m_column_size{column_size_} {
        stats.resize(m_column_size);
    }

    template <typename T, typename TOut>
    void operator()(const T* old_y, const T* old_x1, const T* old_x2, const T* new_y, const T* new_x1, const T* new_x2,
                    TOut* b0, TOut* b1, TOut* b2) {
        for (int i = 0; i < m_column_size; ++i) {
            auto& st = stats[i];
            if (old_y) {
                delete_old(st, old_y[i], old_x1[i], old_x2[i]);
            }
            add_new(st, new_y[i], new_x1[i], new_x2[i]);
            std::tie(b0[i], b1[i], b2[i]) = st.calc();
        }
    }

    template <typename T>
    void delete_old(stat& st, T old_y, T old_x1, T old_x2) {
        if (std::isfinite(old_y) && std::isfinite(old_x1) && std::isfinite(old_x2)) {
            st.sum_x12 -= old_x1 * old_x2;
            st.sum_x1_2 -= old_x1 * old_x1;
            st.sum_x2_2 -= old_x2 * old_x2;
            st.sum_y_2 -= old_y * old_y;
            st.sum_x1y -= old_x1 * old_y;
            st.sum_x2y -= old_x2 * old_y;
            st.sum_x1 -= old_x1;
            st.sum_x2 -= old_x2;
            st.sum_y -= old_y;
            --st.m_valid_count;
        }
    }
    template <typename T>
    void add_new(stat& st, T y, T x1, T x2) {
        if (std::isfinite(y) && std::isfinite(x1) && std::isfinite(x2)) {
            st.sum_x12 += x1 * x2;
            st.sum_x1_2 += x1 * x1;
            st.sum_x2_2 += x2 * x2;
            st.sum_y_2 += y * y;
            st.sum_x1y += x1 * y;
            st.sum_x2y += x2 * y;
            st.sum_x1 += x1;
            st.sum_x2 += x2;
            st.sum_y += y;
            ++st.m_valid_count;
        }
    }

    void init() {
        for (auto& stat : stats) stat.clear();
    }

    template <typename T>
    void full_single(int idx, const T* new_y, const T* new_x1, const T* new_x2) {
        for (int i = 0; i < m_column_size; ++i) {
            auto& st = stats[i];
            add_new(st, new_y[i], new_x1[i], new_x2[i]);
        }
    }

    template <typename TOut>
    void final_result(TOut* b0, TOut* b1, TOut* b2) {
        for (int i = 0; i < m_column_size; ++i) {
            auto& st = stats[i];
            std::tie(b0[i], b1[i], b2[i]) = st.calc();
        }
    }

    void set_row_size(int row) {}
    void set_param(const std::string& key, const std::string& value) {}
};

struct rolling_ema_hl_rb_range {
    struct stat {
        long double total_{0}, total_w{0};
        int m_valid_count{0};
        void clear() {
            total_ = total_w = 0;
            m_valid_count = 0;
        }
        double calc() const {
            if (m_valid_count > 0)
                return total_ / total_w;
            else
                return NAN;
        }
    };
    int m_column_size;
    int m_row_size{0};  // window size
    int count{0};
    double decay_coeff{0};
    std::vector<stat> stats;

    explicit rolling_ema_hl_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    void set_row_size(int row) {
        m_row_size = row;
        decay_coeff = ema_hl2decay(row);
    }
    void set_param(const std::string& key, const std::string& value) {
        if (key == "hl" || key == "half_life") {
            decay_coeff = ema_hl2decay(std::stod(value));
        }
    }

    template <typename T, typename TOut>
    void operator()(const T* old_row, const T* new_row, TOut* output) {
        ++count;
        for (int i = 0; i < m_column_size; ++i) {
            auto& st = stats[i];
            if (old_row) {
                delete_old(st, old_row[i]);
            }
            output[i] = add_new(st, new_row[i]);
        }

        if (count < m_row_size) {
            std::fill(output, output + m_column_size, NAN);
        }
    }

    template <typename T>
    void delete_old(stat& st, T old_value) {
        if (std::isfinite(old_value)) {
            st.total_ -= old_value * std::pow(decay_coeff, m_row_size - 1);
            st.total_w -= std::pow(decay_coeff, m_row_size - 1);
            --st.m_valid_count;
        }
    }

    template <typename T>
    double add_new(stat& st, T new_value) {
        if (std::isfinite(new_value)) {
            st.total_ = st.total_ * decay_coeff + new_value;
            st.total_w = st.total_w * decay_coeff + 1.0;
            ++st.m_valid_count;
        } else {
            st.total_ *= decay_coeff;
            st.total_w *= decay_coeff;
        }
        return st.calc();
    }

    void init() {
        for (auto& stat : stats) stat.clear();
    }

    template <typename T, typename TOut>
    void full_single(int idx, const T* _row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            if (std::isfinite(_row[i])) {
                ++stats[i].m_valid_count;
                stats[i].total_ += _row[i] * pow(decay_coeff, idx);
                stats[i].total_w += pow(decay_coeff, idx);
            }
        }
    }

    template <typename TOut>
    void final_result(TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            output[i] = stats[i].calc();
        }
    }
};

}  // namespace ornate

#endif
