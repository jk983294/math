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

    void set_ins_num(int ins_num) {
        m_column_size = ins_num;
        stats.resize(m_column_size);
    }

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

    double warm_up() {
        double ret = 0.;
        for (auto& stat : stats) ret += stat.m_valid_count;
        return ret / int(stats.size() + 1);
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

    void set_ins_num(int ins_num) {
        m_column_size = ins_num;
        stats.resize(m_column_size);
    }

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

    double warm_up() {
        double ret = 0.;
        for (auto& stat : stats) ret += stat.m_valid_count;
        return ret / int(stats.size() + 1);
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
        int m_valid_count{0}, m_zero_count{0};
        void clear() {
            total_prod = 1;
            m_valid_count = m_zero_count = 0;
        }
        double calc() const {
            if (m_valid_count > 0) {
                if (m_zero_count > 0)
                    return 0;
                else
                    return total_prod;
            } else
                return NAN;
        }
    };
    int m_column_size;
    std::vector<stat> stats;

    explicit rolling_prod_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    void set_ins_num(int ins_num) {
        m_column_size = ins_num;
        stats.resize(m_column_size);
    }

    template <typename T, typename TOut>
    void operator()(const T* old_row, const T* new_row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            auto& stat = stats[i];
            if (old_row) {
                auto old_value = old_row[i];
                if (std::isfinite(old_value)) {
                    if (old_value == 0) {
                        --stat.m_zero_count;
                    } else {
                        stat.total_prod /= old_value;
                    }
                    --stat.m_valid_count;
                }
            }

            auto new_value = new_row[i];
            if (std::isfinite(new_value)) {
                if (new_value == 0) {
                    ++stat.m_zero_count;
                } else {
                    stat.total_prod *= new_value;
                }
                ++stat.m_valid_count;
            }

            output[i] = stat.calc();
        }
    }

    void init() {
        for (auto& stat : stats) stat.clear();
    }

    double warm_up() {
        double ret = 0.;
        for (auto& stat : stats) ret += stat.m_valid_count;
        return ret / int(stats.size() + 1);
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

    void set_ins_num(int ins_num) { m_column_size = ins_num; }

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
    void init() {}
    double warm_up() { return 0.; }
};

struct rolling_delta_rb_range {
    int m_column_size;

    explicit rolling_delta_rb_range(int column_size_) : m_column_size{column_size_} {}

    void set_ins_num(int ins_num) { m_column_size = ins_num; }

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
    void init() {}
    double warm_up() { return 0.; }
};

struct rolling_pct_rb_range {
    int m_column_size;

    explicit rolling_pct_rb_range(int column_size_) : m_column_size{column_size_} {}

    void set_ins_num(int ins_num) { m_column_size = ins_num; }

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
    void init() {}
    double warm_up() { return 0.; }
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
            if (m_valid_count > 2) {
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

    void set_ins_num(int ins_num) {
        m_column_size = ins_num;
        stats.resize(m_column_size);
    }

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

    double warm_up() {
        double ret = 0.;
        for (auto& stat : stats) ret += stat.m_valid_count;
        return ret / int(stats.size() + 1);
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
            if (m_valid_count > 2) {
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

    void set_ins_num(int ins_num) {
        m_column_size = ins_num;
        stats.resize(m_column_size);
    }

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

    double warm_up() {
        double ret = 0.;
        for (auto& stat : stats) ret += stat.m_valid_count;
        return ret / int(stats.size() + 1);
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
            if (m_valid_count > 2 && std::isfinite(latest_val)) {
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

    void set_ins_num(int ins_num) {
        m_column_size = ins_num;
        stats.resize(m_column_size);
    }

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
            if (st.m_valid_count <= 2) return NAN;
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

    double warm_up() {
        double ret = 0.;
        for (auto& stat : stats) ret += stat.m_valid_count;
        return ret / int(stats.size() + 1);
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

    void set_ins_num(int ins_num) {
        m_column_size = ins_num;
        stats.resize(m_column_size);
    }

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

    double warm_up() {
        double ret = 0.;
        for (auto& stat : stats) ret += stat.m_valid_count;
        return ret / int(stats.size() + 1);
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
            if (m_valid_count > 2) {
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

    void set_ins_num(int ins_num) {
        m_column_size = ins_num;
        stats.resize(m_column_size);
    }

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

    double warm_up() {
        double ret = 0.;
        for (auto& stat : stats) ret += stat.m_valid_count;
        return ret / int(stats.size() + 1);
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
            if (m_valid_count > 2) {
                double cov = (sumxy / m_valid_count - sumx * sumy / (m_valid_count * m_valid_count));
                double var1 = (sum_x2 / m_valid_count - sumx * sumx / (m_valid_count * m_valid_count));
                double var2 = (sum_y2 / m_valid_count - sumy * sumy / (m_valid_count * m_valid_count));
                double denominator = sqrt(var1 * var2);
                if (denominator > 0) return cov / denominator;
            }
            return NAN;
        }
    };
    int m_column_size;
    std::vector<stat> stats;

    explicit rolling_corr_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    void set_ins_num(int ins_num) {
        m_column_size = ins_num;
        stats.resize(m_column_size);
    }

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

    double warm_up() {
        double ret = 0.;
        for (auto& stat : stats) ret += stat.m_valid_count;
        return ret / int(stats.size() + 1);
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
            if (m_valid_count > 2) {
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

    void set_ins_num(int ins_num) {
        m_column_size = ins_num;
        stats.resize(m_column_size);
    }

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

    double warm_up() {
        double ret = 0.;
        for (auto& stat : stats) ret += stat.m_valid_count;
        return ret / int(stats.size() + 1);
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
            if (m_valid_count > 2) {
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

    void set_ins_num(int ins_num) {
        m_column_size = ins_num;
        stats.resize(m_column_size);
    }

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

    double warm_up() {
        double ret = 0.;
        for (auto& stat : stats) ret += stat.m_valid_count;
        return ret / int(stats.size() + 1);
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

    void set_ins_num(int ins_num) {
        m_column_size = ins_num;
        stats.resize(m_column_size);
    }

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

    double warm_up() {
        double ret = 0.;
        for (auto& stat : stats) ret += stat.m_valid_count;
        return ret / int(stats.size() + 1);
    }

    template <typename T, typename TOut>
    void full_single(int idx, int window, const T* _row, TOut* output) {
        ++count;
        for (int i = 0; i < m_column_size; ++i) {
            if (std::isfinite(_row[i])) {
                stats[i].m_valid_count += (window - idx);
                stats[i].total += _row[i] * (window - idx);
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
        void clear() {
            front = 0;
            rear = 0;
            seq = -1;
        }
    };
    int capacity{0};
    int m_column_size;
    int m_count{0};
    std::vector<stat> stats;
    std::vector<Cell> _data;
    TCmp cmp;

    explicit rolling_mq_rb_range(int column_size_) : m_column_size{column_size_} {
        stats.resize(m_column_size);
        cmp = TCmp();
    }

    void init() {
        for (auto& st : stats) st.clear();
        m_count = 0;
    }

    double warm_up() {
        double ret = 0.;
        for (auto& stat : stats) ret += stat.seq;
        for (auto& d : _data) ret += d.data;
        return ret / int(stats.size() + _data.size() + 1);
    }

    void set_ins_num(int ins_num) {
        m_column_size = ins_num;
        stats.resize(m_column_size);
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
        if (m_count == 1) {
            return 0.0;
        } else if (m_count < capacity - 1) {
            return double(st.seq - start_cell[st.front].seq) / (m_count - 1);
        } else
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
    using rolling_mq_rb_range<TData, TCmp>::m_count;

    explicit rolling_mq_index_rb_range(int column_size_) : rolling_mq_rb_range<TData, TCmp>(column_size_) {}

    template <typename T, typename TOut>
    void operator()(const T* old_row, const T* new_row, TOut* output) {
        ++m_count;
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
    using rolling_mq_rb_range<TData, TCmp>::m_count;

    explicit rolling_mq_percent_rb_range(int column_size_) : rolling_mq_rb_range<TData, TCmp>(column_size_) {}

    template <typename T, typename TOut>
    void operator()(const T* old_row, const T* new_row, TOut* output) {
        ++m_count;
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
    using rolling_mq_rb_range<TData, TCmp>::m_count;

    explicit rolling_mq_value_rb_range(int column_size_) : rolling_mq_rb_range<TData, TCmp>(column_size_) {}

    template <typename T, typename TOut>
    void operator()(const T* old_row, const T* new_row, TOut* output) {
        ++m_count;
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

    void init() {
        m_count = 0;
        for (auto& st : stats) st.m_valid_count = 0;
    }

    double warm_up() {
        double ret = 0.;
        for (auto& stat : stats) ret += stat.m_valid_count;
        return ret / int(stats.size() + 1);
    }

    void set_ins_num(int ins_num) {
        m_column_size = ins_num;
        stats.resize(m_column_size);
    }

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

    void init() { m_count = 0; }
    double warm_up() {
        long double ret = 0.;
        for (auto& stat : m_container) ret += stat;
        return ret / (m_container.size() + m_count + 1);
    }

    void set_ins_num(int ins_num) { m_column_size = ins_num; }

    void set_row_size(int row) {
        window_size = row;
        m_container.resize(window_size * m_column_size);
    }

    double calc(const T* x, T new_value, int real_size) {
        if (std::isnan(new_value)) {
            return NAN;
        }
        int nlte = 0, neq = 0, nv = 0;

        for (int i = 0; i < real_size; ++i) {
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
        int real_size = std::min(m_count, window_size);
        for (int i = 0; i < m_column_size; ++i) {
            auto* start_data = m_container.data() + i * window_size;
            const T new_val = new_row[i];
            start_data[pos] = new_val;
            output[i] = calc(start_data, new_val, real_size);
        }
    }

    void set_param(const std::string& key, const std::string& value) {}
};

template <typename T>
struct rolling_rank2_count_rb_range {
    int window_size{0}, m_count{0};
    int m_column_size{0};
    std::vector<T> m_container;

    explicit rolling_rank2_count_rb_range(int column_size_) : m_column_size{column_size_} {}

    void set_ins_num(int ins_num) { m_column_size = ins_num; }

    void set_row_size(int row) {
        window_size = row;
        m_container.resize(window_size * m_column_size);
    }

    void init() { m_count = 0; }
    double warm_up() {
        long double ret = 0.;
        for (auto& stat : m_container) ret += stat;
        return ret / (m_container.size() + m_count + window_size + 1);
    }

    double calc(const T* x, T new_value) {
        if (std::isnan(new_value) || m_count < window_size) {
            return NAN;
        }
        int nlte = 0, nv = 0;
        for (int i = 0; i < window_size; ++i) {
            const T val = x[i];
            if (std::isnan(val)) continue;
            if (val <= new_value) nlte++;
            nv++;
        }
        return nv <= 1 ? NAN : (nlte - 0.0) / (nv - 0.0);
    }

    template <typename TOut>
    void operator()(const T* old_cmp_row, const T* old_row, const T* cmp_row, const T* new_row, TOut* output) {
        int pos = m_count % window_size;
        ++m_count;
        for (int i = 0; i < m_column_size; ++i) {
            auto* start_data = m_container.data() + i * window_size;
            start_data[pos] = new_row[i];
            output[i] = calc(start_data, cmp_row[i]);
        }
    }

    void set_param(const std::string& key, const std::string& value) {}
};

struct rolling_regression2_rb_range {
    struct stat {
        long double sum_x{0}, sum_y{0}, sum_x2{0}, sum_xy{0};
        int m_valid_count{0};
        double a{NAN}, b{NAN}, mean_y, res_squared{0}, y_diff_squared{0};
        void clear() {
            a = b = NAN;
            sum_x = sum_y = sum_x2 = sum_xy = 0;
            m_valid_count = 0;
        }
        void calc_coef() {
            a = b = NAN;
            if (m_valid_count > 2) {
                b = (m_valid_count * sum_xy - sum_x * sum_y) / (m_valid_count * sum_x2 - sum_x * sum_x);
                a = (sum_y - b * sum_x) / m_valid_count;
            }
        }
        double calc_fitted(double x) const { return a + b * x; }
        double calc_residual(double x, double y) const { return y - calc_fitted(x); }
        void r2_pre_calc() {
            calc_coef();
            res_squared = y_diff_squared = 0;
            if (m_valid_count > 2) {
                mean_y = sum_y / m_valid_count;
            } else {
                mean_y = NAN;
            }
        }
        void r2_single(double x, double y) {
            if (m_valid_count <= 2) return;
            if (!std::isfinite(x) || !std::isfinite(y)) return;
            double res = calc_residual(x, y);
            res_squared += res * res;
            double diff = y - mean_y;
            y_diff_squared += diff * diff;
        }
        double get_r2() const {
            if (m_valid_count > 2) {
                return 1. - res_squared / y_diff_squared;
            } else {
                return NAN;
            }
        }
    };
    int m_column_size;
    std::vector<stat> stats;

    explicit rolling_regression2_rb_range(int column_size_) : m_column_size{column_size_} {
        stats.resize(m_column_size);
    }

    void set_ins_num(int ins_num) {
        m_column_size = ins_num;
        stats.resize(m_column_size);
    }

    template <typename TOut>
    void get_coefficients(TOut* a, TOut* b) {
        for (int i = 0; i < m_column_size; ++i) {
            stats[i].calc_coef();
            a[i] = stats[i].a;
            b[i] = stats[i].b;
        }
    }

    template <typename T, typename TOut>
    void get_fitted(const T* new_x, TOut* ret) {
        for (int i = 0; i < m_column_size; ++i) {
            stats[i].calc_coef();
            ret[i] = stats[i].calc_fitted(new_x[i]);
        }
    }

    template <typename T, typename TOut>
    void get_residual(const T* new_y, const T* new_x, TOut* ret) {
        for (int i = 0; i < m_column_size; ++i) {
            stats[i].calc_coef();
            ret[i] = stats[i].calc_residual(new_x[i], new_y[i]);
        }
    }

    void r2_pre_calc() {
        for (int i = 0; i < m_column_size; ++i) {
            stats[i].r2_pre_calc();
        }
    }

    template <typename T>
    void r2_single(int idx, const T* y_row, const T* x_row) {
        for (int i = 0; i < m_column_size; ++i) {
            stats[i].r2_single(x_row[i], y_row[i]);
        }
    }

    template <typename TOut>
    void get_r2(TOut* ret) {
        for (int i = 0; i < m_column_size; ++i) {
            ret[i] = stats[i].get_r2();
        }
    }

    template <typename T>
    void operator()(const T* old_y, const T* old_x, const T* new_y, const T* new_x) {
        for (int i = 0; i < m_column_size; ++i) {
            auto& st = stats[i];
            if (old_y) {
                delete_old(st, old_y[i], old_x[i]);
            }
            add_new(st, new_y[i], new_x[i]);
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
    double warm_up() {
        double ret = 0.;
        for (auto& stat : stats) ret += stat.m_valid_count;
        return ret / int(stats.size() + 1);
    }

    template <typename T>
    void full_single(int idx, const T* y_row, const T* x_row) {
        for (int i = 0; i < m_column_size; ++i) {
            auto& st = stats[i];
            add_new(st, y_row[i], x_row[i]);
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
        double b0{NAN}, b1{NAN}, b2{NAN}, mean_y, res_squared{0}, y_diff_squared{0};
        int m_valid_count{0};

        void clear() {
            sum_x1_2 = sum_x2_2 = sum_y_2 = sum_x1 = sum_x2 = sum_y = sum_x12 = sum_x1y = sum_x2y = 0;
            m_valid_count = 0;
        }
        void calc_coef() {
            b0 = b1 = b2 = NAN;
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
        }
        double calc_fitted(double x1, double x2) const { return b0 + b1 * x1 + b2 * x2; }
        double calc_residual(double x1, double x2, double y) const { return y - calc_fitted(x1, x2); }
        void r2_pre_calc() {
            calc_coef();
            res_squared = y_diff_squared = 0;
            mean_y = NAN;
            if (m_valid_count > 1) {
                mean_y = sum_y / m_valid_count;
            }
        }
        void r2_single(double x1, double x2, double y) {
            if (m_valid_count <= 1) return;
            if (!std::isfinite(x1) || !std::isfinite(x2) || !std::isfinite(y)) return;
            double res = calc_residual(x1, x2, y);
            res_squared += res * res;
            double diff = y - mean_y;
            y_diff_squared += diff * diff;
        }
        double get_r2() const {
            if (m_valid_count > 1) {
                return 1. - res_squared / y_diff_squared;
            } else {
                return NAN;
            }
        }
    };
    int m_column_size;
    std::vector<stat> stats;

    explicit rolling_regression3_rb_range(int column_size_) : m_column_size{column_size_} {
        stats.resize(m_column_size);
    }

    void set_ins_num(int ins_num) {
        m_column_size = ins_num;
        stats.resize(m_column_size);
    }

    template <typename TOut>
    void get_coefficients(TOut* b0, TOut* b1, TOut* b2) {
        for (int i = 0; i < m_column_size; ++i) {
            stats[i].calc_coef();
            b0[i] = stats[i].b0;
            b1[i] = stats[i].b1;
            b2[i] = stats[i].b2;
        }
    }

    template <typename T, typename TOut>
    void get_fitted(const T* new_x1, const T* new_x2, TOut* ret) {
        for (int i = 0; i < m_column_size; ++i) {
            stats[i].calc_coef();
            ret[i] = stats[i].calc_fitted(new_x1[i], new_x2[i]);
        }
    }

    template <typename T, typename TOut>
    void get_residual(const T* new_y, const T* new_x1, const T* new_x2, TOut* ret) {
        for (int i = 0; i < m_column_size; ++i) {
            stats[i].calc_coef();
            ret[i] = stats[i].calc_residual(new_x1[i], new_x2[i], new_y[i]);
        }
    }

    void r2_pre_calc() {
        for (int i = 0; i < m_column_size; ++i) {
            stats[i].r2_pre_calc();
        }
    }

    template <typename T>
    void r2_single(int idx, const T* y_row, const T* x1_row, const T* x2_row) {
        for (int i = 0; i < m_column_size; ++i) {
            stats[i].r2_single(x1_row[i], x2_row[i], y_row[i]);
        }
    }

    template <typename TOut>
    void get_r2(TOut* ret) {
        for (int i = 0; i < m_column_size; ++i) {
            ret[i] = stats[i].get_r2();
        }
    }

    template <typename T>
    void operator()(const T* old_y, const T* old_x1, const T* old_x2, const T* new_y, const T* new_x1,
                    const T* new_x2) {
        for (int i = 0; i < m_column_size; ++i) {
            auto& st = stats[i];
            if (old_y) {
                delete_old(st, old_y[i], old_x1[i], old_x2[i]);
            }
            add_new(st, new_y[i], new_x1[i], new_x2[i]);
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
    double warm_up() {
        double ret = 0.;
        for (auto& stat : stats) ret += stat.m_valid_count;
        return ret / int(stats.size() + 1);
    }

    template <typename T>
    void full_single(int idx, const T* new_y, const T* new_x1, const T* new_x2) {
        for (int i = 0; i < m_column_size; ++i) {
            auto& st = stats[i];
            add_new(st, new_y[i], new_x1[i], new_x2[i]);
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
    double remove_coeff{0};
    std::vector<double> cached_coeff;
    std::vector<stat> stats;

    explicit rolling_ema_hl_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    void set_ins_num(int ins_num) {
        m_column_size = ins_num;
        stats.resize(m_column_size);
    }

    void set_row_size(int row) {
        m_row_size = row;
        decay_coeff = ema_hl2decay(row);
        calc_cached_coeff();
    }
    void set_param(const std::string& key, const std::string& value) {
        if (key == "hl" || key == "half_life" || key == "arg1") {
            decay_coeff = ema_hl2decay(std::stod(value));
            calc_cached_coeff();
        }
    }

    void calc_cached_coeff() {
        remove_coeff = std::pow(decay_coeff, m_row_size - 1);
        cached_coeff.resize(m_row_size);
        for (int i = 0; i < m_row_size; ++i) {
            cached_coeff[i] = std::pow(decay_coeff, i);
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
    }

    template <typename T>
    void delete_old(stat& st, T old_value) {
        if (std::isfinite(old_value)) {
            st.total_ -= old_value * remove_coeff;
            st.total_w -= remove_coeff;
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
    double warm_up() {
        double ret = 0.;
        for (auto& stat : stats) ret += stat.m_valid_count;
        return ret / int(stats.size() + 1);
    }

    template <typename T, typename TOut>
    void full_single(int idx, const T* _row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            if (std::isfinite(_row[i])) {
                ++stats[i].m_valid_count;
                stats[i].total_ += _row[i] * cached_coeff[idx];
                stats[i].total_w += cached_coeff[idx];
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

struct rolling_ema_hl2_rb_range {
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
    double remove_coeff{0};
    std::vector<double> cached_coeff;
    double decay_coeff2{0};
    double remove_coeff2{0};
    std::vector<double> cached_coeff2;
    std::vector<stat> stats;

    explicit rolling_ema_hl2_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    void set_ins_num(int ins_num) {
        m_column_size = ins_num;
        stats.resize(m_column_size);
    }

    void set_row_size(int row) { m_row_size = row; }
    void set_param(const std::string& key, const std::string& value) {
        if (key == "arg1") {
            decay_coeff = ema_hl2decay(std::stod(value));
            remove_coeff = std::pow(decay_coeff, m_row_size - 1);
            cached_coeff.resize(m_row_size);
            for (int i = 0; i < m_row_size; ++i) {
                cached_coeff[i] = std::pow(decay_coeff, i);
            }
        } else if (key == "arg2") {
            decay_coeff2 = ema_hl2decay(std::stod(value));
            remove_coeff2 = std::pow(decay_coeff2, m_row_size - 1);
            cached_coeff2.resize(m_row_size);
            for (int i = 0; i < m_row_size; ++i) {
                cached_coeff2[i] = std::pow(decay_coeff2, i);
            }
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
    }

    template <typename T>
    void delete_old(stat& st, T old_value) {
        if (std::isfinite(old_value)) {
            st.total_ -= old_value * remove_coeff;
            st.total_w -= remove_coeff2;
            --st.m_valid_count;
        }
    }

    template <typename T>
    double add_new(stat& st, T new_value) {
        if (std::isfinite(new_value)) {
            st.total_ = st.total_ * decay_coeff + new_value;
            st.total_w = st.total_w * decay_coeff2 + 1.0;
            ++st.m_valid_count;
        } else {
            st.total_ *= decay_coeff;
            st.total_w *= decay_coeff2;
        }
        return st.calc();
    }

    void init() {
        for (auto& stat : stats) stat.clear();
    }
    double warm_up() {
        double ret = 0.;
        for (auto& stat : stats) ret += stat.m_valid_count;
        return ret / int(stats.size() + 1);
    }

    template <typename T, typename TOut>
    void full_single(int idx, const T* _row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            if (std::isfinite(_row[i])) {
                ++stats[i].m_valid_count;
                stats[i].total_ += _row[i] * cached_coeff[idx];
                stats[i].total_w += cached_coeff2[idx];
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

struct rolling_ols2_rb_range {
    struct stat {
        long double sum_x2{0}, sum_xy{0};
        int m_valid_count{0};
        double b{NAN}, res_squared{0}, sum_y2{0};
        void clear() {
            sum_x2 = sum_xy = 0;
            m_valid_count = 0;
            b = NAN;
        }

        void calc_coef() {
            b = NAN;
            if (m_valid_count > 2) {
                b = sum_xy / sum_x2;
            }
        }
        double calc_fitted(double x) const { return b * x; }
        double calc_residual(double x, double y) const { return y - calc_fitted(x); }
        void r2_pre_calc() {
            calc_coef();
            res_squared = sum_y2 = 0;
        }
        void r2_single(double x, double y) {
            if (m_valid_count < 3) return;
            if (!std::isfinite(x) || !std::isfinite(y)) return;
            double res = calc_residual(x, y);
            res_squared += res * res;
            sum_y2 += y * y;
        }
        double get_r2() const {
            if (m_valid_count > 2) {
                return 1. - res_squared / sum_y2;
            } else {
                return NAN;
            }
        }
    };
    int m_column_size;
    std::vector<stat> stats;

    explicit rolling_ols2_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    void set_ins_num(int ins_num) {
        m_column_size = ins_num;
        stats.resize(m_column_size);
    }

    template <typename TOut>
    void get_coefficients(TOut* b) {
        for (int i = 0; i < m_column_size; ++i) {
            stats[i].calc_coef();
            b[i] = stats[i].b;
        }
    }

    template <typename T, typename TOut>
    void get_fitted(const T* new_x, TOut* ret) {
        for (int i = 0; i < m_column_size; ++i) {
            stats[i].calc_coef();
            ret[i] = stats[i].calc_fitted(new_x[i]);
        }
    }

    template <typename T, typename TOut>
    void get_residual(const T* new_y, const T* new_x, TOut* ret) {
        for (int i = 0; i < m_column_size; ++i) {
            stats[i].calc_coef();
            ret[i] = stats[i].calc_residual(new_x[i], new_y[i]);
        }
    }

    void r2_pre_calc() {
        for (int i = 0; i < m_column_size; ++i) {
            stats[i].r2_pre_calc();
        }
    }

    template <typename T>
    void r2_single(int idx, const T* y_row, const T* x_row) {
        for (int i = 0; i < m_column_size; ++i) {
            stats[i].r2_single(x_row[i], y_row[i]);
        }
    }

    template <typename TOut>
    void get_r2(TOut* ret) {
        for (int i = 0; i < m_column_size; ++i) {
            ret[i] = stats[i].get_r2();
        }
    }

    template <typename T>
    void operator()(const T* old_y, const T* old_x, const T* new_y, const T* new_x) {
        for (int i = 0; i < m_column_size; ++i) {
            auto& st = stats[i];
            if (old_y) {
                delete_old(st, old_y[i], old_x[i]);
            }
            add_new(st, new_y[i], new_x[i]);
        }
    }

    template <typename T>
    void delete_old(stat& st, T old_y, T old_x) {
        if (std::isfinite(old_y) && std::isfinite(old_x)) {
            st.sum_x2 -= old_x * old_x;
            st.sum_xy -= old_y * old_x;
            --st.m_valid_count;
        }
    }
    template <typename T>
    void add_new(stat& st, T y, T x) {
        if (std::isfinite(y) && std::isfinite(x)) {
            st.sum_x2 += x * x;
            st.sum_xy += y * x;
            ++st.m_valid_count;
        }
    }

    void init() {
        for (auto& stat : stats) stat.clear();
    }
    double warm_up() {
        double ret = 0.;
        for (auto& stat : stats) ret += stat.m_valid_count;
        return ret / int(stats.size() + 1);
    }

    template <typename T>
    void full_single(int idx, const T* y_row, const T* x_row) {
        for (int i = 0; i < m_column_size; ++i) {
            auto& st = stats[i];
            add_new(st, y_row[i], x_row[i]);
        }
    }

    void set_row_size(int row) {}
    void set_param(const std::string& key, const std::string& value) {}
};

struct rolling_ols3_rb_range {
    struct stat {
        long double sum_x12 = 0, sum_x1_2 = 0, sum_x2_2 = 0, sum_x1y = 0, sum_x2y = 0;
        int m_valid_count{0};
        double b1{NAN}, b2{NAN}, res_squared{0}, y_squared{0};
        void clear() {
            sum_x1_2 = sum_x2_2 = sum_x12 = sum_x1y = sum_x2y = 0;
            m_valid_count = 0;
        }
        void calc_coef() {
            b1 = b2 = NAN;
            if (m_valid_count >= 2) {
                long double denominator = sum_x12 * sum_x12 - sum_x1_2 * sum_x2_2;
                b1 = (sum_x2y * sum_x12 - sum_x1y * sum_x2_2) / denominator;
                b2 = (sum_x1y * sum_x12 - sum_x2y * sum_x1_2) / denominator;
            }
        }
        double calc_fitted(double x1, double x2) const { return b1 * x1 + b2 * x2; }
        double calc_residual(double x1, double x2, double y) const { return y - calc_fitted(x1, x2); }
        void r2_pre_calc() {
            calc_coef();
            res_squared = y_squared = 0;
        }
        void r2_single(double x1, double x2, double y) {
            if (m_valid_count <= 1) return;
            if (!std::isfinite(x1) || !std::isfinite(x2) || !std::isfinite(y)) return;
            double res = calc_residual(x1, x2, y);
            res_squared += res * res;
            y_squared += y * y;
        }
        double get_r2() const {
            if (m_valid_count >= 2) {
                return 1. - res_squared / y_squared;
            } else {
                return NAN;
            }
        }
    };
    int m_column_size;
    std::vector<stat> stats;

    explicit rolling_ols3_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    void set_ins_num(int ins_num) {
        m_column_size = ins_num;
        stats.resize(m_column_size);
    }

    template <typename TOut>
    void get_coefficients(TOut* b1, TOut* b2) {
        for (int i = 0; i < m_column_size; ++i) {
            stats[i].calc_coef();
            b1[i] = stats[i].b1;
            b2[i] = stats[i].b2;
        }
    }

    template <typename T, typename TOut>
    void get_fitted(const T* new_x1, const T* new_x2, TOut* ret) {
        for (int i = 0; i < m_column_size; ++i) {
            stats[i].calc_coef();
            ret[i] = stats[i].calc_fitted(new_x1[i], new_x2[i]);
        }
    }

    template <typename T, typename TOut>
    void get_residual(const T* new_y, const T* new_x1, const T* new_x2, TOut* ret) {
        for (int i = 0; i < m_column_size; ++i) {
            stats[i].calc_coef();
            ret[i] = stats[i].calc_residual(new_x1[i], new_x2[i], new_y[i]);
        }
    }

    void r2_pre_calc() {
        for (int i = 0; i < m_column_size; ++i) {
            stats[i].r2_pre_calc();
        }
    }

    template <typename T>
    void r2_single(int idx, const T* y_row, const T* x1_row, const T* x2_row) {
        for (int i = 0; i < m_column_size; ++i) {
            stats[i].r2_single(x1_row[i], x2_row[i], y_row[i]);
        }
    }

    template <typename TOut>
    void get_r2(TOut* ret) {
        for (int i = 0; i < m_column_size; ++i) {
            ret[i] = stats[i].get_r2();
        }
    }

    template <typename T>
    void operator()(const T* old_y, const T* old_x1, const T* old_x2, const T* new_y, const T* new_x1,
                    const T* new_x2) {
        for (int i = 0; i < m_column_size; ++i) {
            auto& st = stats[i];
            if (old_y) {
                delete_old(st, old_y[i], old_x1[i], old_x2[i]);
            }
            add_new(st, new_y[i], new_x1[i], new_x2[i]);
        }
    }

    template <typename T>
    void delete_old(stat& st, T old_y, T old_x1, T old_x2) {
        if (std::isfinite(old_y) && std::isfinite(old_x1) && std::isfinite(old_x2)) {
            st.sum_x12 -= old_x1 * old_x2;
            st.sum_x1_2 -= old_x1 * old_x1;
            st.sum_x2_2 -= old_x2 * old_x2;
            st.sum_x1y -= old_x1 * old_y;
            st.sum_x2y -= old_x2 * old_y;
            --st.m_valid_count;
        }
    }
    template <typename T>
    void add_new(stat& st, T y, T x1, T x2) {
        if (std::isfinite(y) && std::isfinite(x1) && std::isfinite(x2)) {
            st.sum_x12 += x1 * x2;
            st.sum_x1_2 += x1 * x1;
            st.sum_x2_2 += x2 * x2;
            st.sum_x1y += x1 * y;
            st.sum_x2y += x2 * y;
            ++st.m_valid_count;
        }
    }

    void init() {
        for (auto& stat : stats) stat.clear();
    }
    double warm_up() {
        double ret = 0.;
        for (auto& stat : stats) ret += stat.m_valid_count;
        return ret / int(stats.size() + 1);
    }

    template <typename T>
    void full_single(int idx, const T* new_y, const T* new_x1, const T* new_x2) {
        for (int i = 0; i < m_column_size; ++i) {
            auto& st = stats[i];
            add_new(st, new_y[i], new_x1[i], new_x2[i]);
        }
    }

    void set_row_size(int row) {}
    void set_param(const std::string& key, const std::string& value) {}
};

/**
 * 因为 sum_x2 不能 rolling, 所以没有 rolling 算法
 * 如果 input 没有 nan,那么可以rolling, 当有nan数据的时候, sigma(i*i)没法rolling
 */
struct slope_no_intercept_rb_range {
    struct stat {
        long double sum_x2{0}, sum_xy{0};
        int m_valid_count{0};
        void clear() {
            sum_x2 = sum_xy = 0;
            m_valid_count = 0;
        }
        double calc() const {
            if (m_valid_count > 2) {
                return sum_xy / sum_x2;
            }
            return NAN;
        }
    };
    int m_column_size;
    std::vector<stat> stats;
    std::vector<double> x_squared;
    int m_row_size{0};  // window size

    explicit slope_no_intercept_rb_range(int column_size_) : m_column_size{column_size_} {
        stats.resize(m_column_size);
    }

    void set_ins_num(int ins_num) {
        m_column_size = ins_num;
        stats.resize(m_column_size);
    }

    template <typename T>
    void add_new(stat& st, T y, int idx_) {
        if (std::isfinite(y)) {
            st.sum_x2 += x_squared[idx_];
            st.sum_xy += y * idx_;
            ++st.m_valid_count;
        }
    }

    void init() {
        for (auto& stat : stats) stat.clear();
    }
    double warm_up() {
        double ret = 0.;
        for (auto& stat : stats) ret += stat.m_valid_count;
        return ret / int(stats.size() + 1);
    }

    template <typename T, typename TOut>
    void full_single(int idx, int window, const T* y_row, TOut* output) {
        int x_idx = window - 1 - idx;
        for (int i = 0; i < m_column_size; ++i) {
            auto& st = stats[i];
            add_new(st, y_row[i], x_idx);
        }
    }

    template <typename TOut>
    void final_result(TOut* b) {
        for (int i = 0; i < m_column_size; ++i) {
            auto& st = stats[i];
            b[i] = st.calc();
        }
    }

    void set_row_size(int row) {
        m_row_size = row;
        x_squared.resize(m_row_size);
        for (int i = 0; i < m_row_size; ++i) {
            x_squared[i] = i * i;
        }
    }
    void set_param(const std::string& key, const std::string& value) {}
};

struct slope_rb_range {
    struct stat {
        long double sum_x2{0}, sum_xy{0}, sum_x{0}, sum_y{0};
        int m_valid_count{0};
        void clear() {
            sum_x2 = sum_xy = sum_x = sum_y = 0;
            m_valid_count = 0;
        }
        double calc() const {
            if (m_valid_count > 2) {
                long double cov_xy = sum_xy * m_valid_count - sum_x * sum_y;
                long double var_x = sum_x2 * m_valid_count - sum_x * sum_x;
                return var_x > 0 ? cov_xy / var_x : NAN;
            }
            return NAN;
        }
    };
    int m_column_size;
    std::vector<stat> stats;
    std::vector<double> x_squared;
    int m_row_size{0};  // window size

    explicit slope_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    void set_ins_num(int ins_num) {
        m_column_size = ins_num;
        stats.resize(m_column_size);
    }

    template <typename T>
    void add_new(stat& st, T y, int idx_) {
        if (std::isfinite(y)) {
            st.sum_x2 += x_squared[idx_];
            st.sum_xy += y * idx_;
            st.sum_x += idx_;
            st.sum_y += y;
            ++st.m_valid_count;
        }
    }

    void init() {
        for (auto& stat : stats) stat.clear();
    }
    double warm_up() {
        double ret = 0.;
        for (auto& stat : stats) ret += stat.m_valid_count;
        return ret / int(stats.size() + 1);
    }

    template <typename T, typename TOut>
    void full_single(int idx, int window, const T* y_row, TOut* output) {
        int x_idx = window - 1 - idx;
        for (int i = 0; i < m_column_size; ++i) {
            auto& st = stats[i];
            add_new(st, y_row[i], x_idx);
        }
    }

    template <typename TOut>
    void final_result(TOut* b) {
        for (int i = 0; i < m_column_size; ++i) {
            auto& st = stats[i];
            b[i] = st.calc();
        }
    }

    void set_row_size(int row) {
        m_row_size = row;
        x_squared.resize(m_row_size);
        for (int i = 0; i < m_row_size; ++i) {
            x_squared[i] = i * i;
        }
    }
    void set_param(const std::string& key, const std::string& value) {}
};

struct rolling_sharpe_rb_range {
    struct stat {
        long double total_sum{0}, total_square_sum{0};
        int m_valid_count{0};
        void clear() {
            total_sum = 0;
            total_square_sum = 0;
            m_valid_count = 0;
        }
        double calc() const {
            if (m_valid_count <= 2) return NAN;
            long double mean = total_sum / m_valid_count;
            long double sd = sqrt((total_square_sum - mean * mean * m_valid_count) / (m_valid_count - 1.0));
            return sd > 1e-16 ? mean / sd : NAN;
        }
    };
    int m_column_size;
    std::vector<stat> stats;

    explicit rolling_sharpe_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    void set_ins_num(int ins_num) {
        m_column_size = ins_num;
        stats.resize(m_column_size);
    }

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
        return st.calc();
    }

    void init() {
        for (auto& stat : stats) stat.clear();
    }
    double warm_up() {
        double ret = 0.;
        for (auto& stat : stats) ret += stat.m_valid_count;
        return ret / int(stats.size() + 1);
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
            output[i] = stats[i].calc();
        }
    }

    void set_row_size(int row) {}
    void set_param(const std::string& key, const std::string& value) {}
};

struct rolling_scale_rb_range {
    struct stat {
        long double total_abs_sum{0};
        int m_valid_count{0};
        void clear() {
            total_abs_sum = 0;
            m_valid_count = 0;
        }
        double calc(double latest_val) const {
            if (std::isfinite(latest_val)) {
                double abs_mean = total_abs_sum / m_valid_count;
                return abs_mean > 1e-6 ? (latest_val / abs_mean) : NAN;
            }
            return NAN;
        }
    };
    int m_column_size;
    std::vector<stat> stats;

    explicit rolling_scale_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    void set_ins_num(int ins_num) {
        m_column_size = ins_num;
        stats.resize(m_column_size);
    }

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
            st.total_abs_sum -= std::abs(old_value);
            --st.m_valid_count;
        }
    }

    template <typename T>
    double add_new(stat& st, T new_value) {
        if (std::isfinite(new_value)) {
            st.total_abs_sum += std::abs(new_value);
            ++st.m_valid_count;
            double abs_mean = st.total_abs_sum / st.m_valid_count;
            return abs_mean > 1e-6 ? (new_value / abs_mean) : NAN;
        } else {
            return NAN;
        }
    }

    void init() {
        for (auto& stat : stats) stat.clear();
    }
    double warm_up() {
        double ret = 0.;
        for (auto& stat : stats) ret += stat.m_valid_count;
        return ret / (stats.size() + 1);
    }

    template <typename T, typename TOut>
    void full_single(int idx, const T* _row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            if (std::isfinite(_row[i])) {
                ++stats[i].m_valid_count;
                stats[i].total_abs_sum += std::abs(_row[i]);
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

struct rolling_rsharpe_rb_range {  // reverse sharpe = sd / mean
    struct stat {
        long double total_sum{0}, total_square_sum{0};
        int m_valid_count{0};
        void clear() {
            total_sum = 0;
            total_square_sum = 0;
            m_valid_count = 0;
        }
        double calc() const {
            if (m_valid_count <= 2) return NAN;
            double mean = total_sum / m_valid_count;
            double sd = sqrt((total_square_sum - mean * mean * m_valid_count) / (m_valid_count - 1.0));
            return sd > 0 ? sd / mean : NAN;
        }
    };
    int m_column_size;
    std::vector<stat> stats;

    explicit rolling_rsharpe_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    void set_ins_num(int ins_num) {
        m_column_size = ins_num;
        stats.resize(m_column_size);
    }

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
        return st.calc();
    }

    void init() {
        for (auto& stat : stats) stat.clear();
    }
    double warm_up() {
        double ret = 0.;
        for (auto& stat : stats) ret += stat.m_valid_count;
        return ret / int(stats.size() + 1);
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
            output[i] = stats[i].calc();
        }
    }

    void set_row_size(int row) {}
    void set_param(const std::string& key, const std::string& value) {}
};

struct rolling_dcor_rb_range {
    struct stat {
        long double sumxy{0}, sum_x2{0}, sum_y2{0};
        int m_valid_count{0};
        void clear() {
            sumxy = sum_x2 = sum_y2 = 0;
            m_valid_count = 0;
        }
        double calc() const {
            if (m_valid_count > 2) {
                long double d = sum_x2 * sum_y2;
                return d > 1e-16 ? sumxy / std::sqrt(d) : NAN;
            }
            return NAN;
        }
    };
    int m_column_size;
    std::vector<stat> stats;

    explicit rolling_dcor_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    void set_ins_num(int ins_num) {
        m_column_size = ins_num;
        stats.resize(m_column_size);
    }

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
            st.sum_x2 -= old_value0 * old_value0;
            st.sum_y2 -= old_value1 * old_value1;
            --st.m_valid_count;
        }
    }

    template <typename T>
    double add_new(stat& st, T data0, T data1) {
        if (std::isfinite(data0) && std::isfinite(data1)) {
            st.sumxy += data0 * data1;
            st.sum_x2 += data0 * data0;
            st.sum_y2 += data1 * data1;
            ++st.m_valid_count;
        }
        return st.calc();
    }

    void init() {
        for (auto& stat : stats) stat.clear();
    }
    double warm_up() {
        double ret = 0.;
        for (auto& stat : stats) ret += stat.m_valid_count;
        return ret / int(stats.size() + 1);
    }

    template <typename T, typename TOut>
    void full_single(int idx, const T* x_row, const T* y_row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            if (std::isfinite(x_row[i]) && std::isfinite(y_row[i])) {
                ++stats[i].m_valid_count;
                stats[i].sumxy += x_row[i] * y_row[i];
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

struct ts_cross_rb_range {
    int m_column_size;

    explicit ts_cross_rb_range(int column_size_) : m_column_size{column_size_} {}

    void set_ins_num(int ins_num) { m_column_size = ins_num; }

    template <typename T, typename TOut>
    void operator()(const T* _row_x0, const T* _row_y0, const T* _row_x1, const T* _row_y1, TOut* output) {
        if (_row_x0 == nullptr || _row_y0 == nullptr) {
            std::fill(output, output + m_column_size, NAN);
            return;
        }
        for (int i = 0; i < m_column_size; ++i) {
            double x0 = _row_x0[i];
            double y0 = _row_y0[i];
            double x1 = _row_x1[i];
            double y1 = _row_y1[i];
            if (std::isnan(x0) || std::isnan(y0) || std::isnan(x1) || std::isnan(y1)) {
                output[i] = NAN;
            } else if (x0 < y0 && x1 > y1) {
                output[i] = 1;
            } else if (x0 > y0 && x1 < y1) {
                output[i] = -1;
            } else {
                output[i] = 0;
            }
        }
    }

    void set_row_size(int row) {}
    void set_param(const std::string& key, const std::string& value) {}
    double warm_up() { return 0.; }
};

/**
 * 没有rolling算法，因为要按当前值来看是否 count++, 下一次，依照的当前值变化了
 * 返回从后往前和当前值同号的数字个数，一旦发现不符合的就结束累计
 */
struct rolling_backward_cpn_rb_range {
    struct stat {
        int m_count{0};
        bool m_is_finish{false};
        void clear() {
            m_is_finish = false;
            m_count = 0;
        }
        double calc() const { return m_count; }
    };
    int m_column_size;
    std::vector<stat> stats;
    int sign{0};  // 0: all, 1: + only, -1: - only

    explicit rolling_backward_cpn_rb_range(int column_size_) : m_column_size{column_size_} {
        stats.resize(m_column_size);
    }

    void set_ins_num(int ins_num) {
        m_column_size = ins_num;
        stats.resize(m_column_size);
    }

    bool should_account(double value, double latest_value) const {
        return (sign != 0 && is_same_sign(value, sign)) || is_same_sign(value, latest_value);
    }

    void init() {
        for (auto& stat : stats) stat.clear();
    }
    double warm_up() {
        double ret = 0.;
        for (auto& stat : stats) ret += stat.m_count;
        return ret / int(stats.size() + 1);
    }

    template <typename T, typename TOut>
    void full_single(int idx, const T* _row, TOut* output) {
        if (idx == 0) {
            detail::_data_copy2vector(_row, output, m_column_size);
        }

        for (int i = 0; i < m_column_size; ++i) {
            if (std::isfinite(_row[i]) && std::isfinite(output[i])) {
                if (!stats[i].m_is_finish && should_account(_row[i], output[i])) {
                    ++stats[i].m_count;
                } else {
                    stats[i].m_is_finish = true;
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

/**
 * 返回window内 x[i] 和 x[i+lag] 同号的比例
 */
struct rolling_ts_acp_rb_range {
    struct stat {
        int m_count{0}, m_pos_count{0};
        void clear() {
            m_count = 0;
            m_pos_count = 0;
        }
        double calc() const { return m_pos_count * 1.0 / m_count; }
    };
    int m_column_size;
    std::vector<stat> stats;

    explicit rolling_ts_acp_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    void set_ins_num(int ins_num) {
        m_column_size = ins_num;
        stats.resize(m_column_size);
    }

    void init() {
        for (auto& stat : stats) stat.clear();
    }
    double warm_up() {
        double ret = 0.;
        for (auto& stat : stats) ret += stat.m_count;
        return ret / int(stats.size() + 1);
    }

    template <typename T>
    void full_single(int idx, const T* _row, const T* _lag_row) {
        for (int i = 0; i < m_column_size; ++i) {
            if (std::isfinite(_row[i]) && std::isfinite(_lag_row[i])) {
                ++stats[i].m_count;
                if (is_same_sign(_row[i], _lag_row[i])) {
                    ++stats[i].m_pos_count;
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
    void set_param(const std::string& key, const std::string& value) {}
};

struct rolling_beta_rb_range {
    struct stat {
        long double sumxy{0}, sum_x2{0}, sum_x{0}, sum_y{0};
        int m_valid_count{0};
        void clear() {
            sumxy = sum_x2 = sum_x = sum_y = 0;
            m_valid_count = 0;
        }
        double calc() const {
            if (m_valid_count > 2) {
                long double d = m_valid_count * sum_x2 - sum_x * sum_y;
                return d > 1e-16 ? (m_valid_count * sumxy - sum_x * sum_y) / d : NAN;
            }
            return NAN;
        }
    };
    int m_column_size;
    std::vector<stat> stats;

    explicit rolling_beta_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    void set_ins_num(int ins_num) {
        m_column_size = ins_num;
        stats.resize(m_column_size);
    }

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
            st.sum_x2 -= old_value0 * old_value0;
            st.sum_x -= old_value0;
            st.sum_y -= old_value1;
            --st.m_valid_count;
        }
    }

    template <typename T>
    double add_new(stat& st, T data0, T data1) {
        if (std::isfinite(data0) && std::isfinite(data1)) {
            st.sumxy += data0 * data1;
            st.sum_x2 += data0 * data0;
            st.sum_x += data0;
            st.sum_y += data1;
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
                stats[i].sum_x2 += x_row[i] * x_row[i];
                stats[i].sum_x += x_row[i];
                stats[i].sum_y += y_row[i];
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

struct CondMean {
    double res = 0;
    int n = 0;
    void operator()(bool c, double val) {
        if (c && std::isfinite(val)) {
            res += val;
            ++n;
        }
    }
    double calc(double fill) { return n > 0 ? res / n : fill; }
};

struct CondSd {
    double sx = 0, sxx = 0;
    int n = 0;
    void operator()(bool c, double val) {
        if (c && std::isfinite(val)) {
            sx += val;
            sxx += val * val;
            ++n;
        }
    }
    double calc(double fill) { return n >= 2 ? sqrt((sxx / n - sx * sx / (n * n)) * n / (n - 1)) : fill; }
};

struct CondMax {
    double res = std::numeric_limits<double>::lowest();
    int n = 0;
    void operator()(bool c, double val) {
        if (c && val > res) {
            res = val;
            ++n;
        }
    }
    double calc(double fill) { return n > 0 ? res : fill; }
};

struct CondMin {
    double res = std::numeric_limits<double>::max();
    int n = 0;
    void operator()(bool c, double val) {
        if (c && val < res) {
            res = val;
            ++n;
        }
    }
    double calc(double fill) { return n > 0 ? res : fill; }
};

template <typename T, typename TComp, typename TCond>
struct rolling_cond_stat_rb_range {
public:
    int window_size{0}, m_count{0};
    int m_column_size;
    double q = 0.5, fill = NAN;
    int method = 1;
    int least = 3;
    bool partial = false;
    std::vector<double> m_cond_data_;
    std::vector<double> m_val_data_;
    std::vector<double> internal_;
    TComp g;

    explicit rolling_cond_stat_rb_range(int column_size_) : m_column_size{column_size_} {}

    double mean(const double* x, int start, int end) {
        double res = 0, n = 0;
        for (int i = start; i < end; ++i) {
            int idx = i % window_size;
            if (std::isfinite(x[idx])) {
                res += x[idx];
                n++;
            }
        }
        return n > 0 ? res / n : fill;
    }

    double quantile(const double* x, int start, int end) {
        double res;
        std::size_t ny = 0;
        for (int i = start; i < end; ++i) {
            int idx = i % window_size;
            if (std::isfinite(x[idx])) {
                internal_[ny++] = x[idx];
            }
        }

        if (ny == 0) {
            return fill;
        }

        double idx = (ny - 1) * q;
        double idx_lb = std::floor(idx);
        double idx_ub = std::ceil(idx);
        if (idx_lb == idx_ub) {
            std::nth_element(internal_.begin(), internal_.begin() + idx, internal_.begin() + ny);
            res = internal_[idx];
        } else {
            std::nth_element(internal_.begin(), internal_.begin() + idx_ub, internal_.begin() + ny);
            std::nth_element(internal_.begin(), internal_.begin() + idx_lb, internal_.begin() + idx_ub);
            res = internal_[idx_lb] * (idx_ub - idx) + internal_[idx_ub] * (idx - idx_lb);
        }
        return res;
    }

    double max(const double* x, int start, int end) {
        double res = std::numeric_limits<double>::lowest();
        int n = 0;
        for (int i = start; i < end; ++i) {
            int idx = i % window_size;
            if (std::isfinite(x[idx]) && x[idx] > res) {
                res = x[idx];
                ++n;
            }
        }
        return n > 0 ? res : fill;
    }

    double min(const double* x, int start, int end) {
        double res = std::numeric_limits<double>::max();
        int n = 0;
        for (int i = start; i < end; ++i) {
            int idx = i % window_size;
            if (std::isfinite(x[idx]) && x[idx] < res) {
                res = x[idx];
                ++n;
            }
        }
        return n > 0 ? res : fill;
    }

    void init() { m_count = 0; }

    double warm_up() { return 0; }

    void set_ins_num(int ins_num) { m_column_size = ins_num; }

    void set_row_size(int row) {
        window_size = row;
        m_cond_data_.resize(row * m_column_size, NAN);
        m_val_data_.resize(row * m_column_size, NAN);
        internal_.resize(window_size);
    }

    double handle(T cond_value, T val_value, double* start_cond_data, double* start_val_data) {
        start_cond_data[(m_count - 1) % window_size] = cond_value;
        start_val_data[(m_count - 1) % window_size] = val_value;
        if (partial && least > 0 && m_count < least) return NAN;
        if (!partial && m_count < window_size) return NAN;
        double xx = cond_value;
        int end = m_count;
        int start = end - window_size;
        if (start < 0) start = 0;

        if (method == 1) {
            xx = mean(start_cond_data, start, end);
        } else if (method == 2) {
            xx = quantile(start_cond_data, start, end);
        } else if (method == 3) {
            xx = max(start_cond_data, start, end);
        } else if (method == 4) {
            xx = min(start_cond_data, start, end);
        }

        TCond stat;
        for (int i = start; i < end; ++i) {
            int idx = i % window_size;
            stat(g(start_cond_data[idx], xx), start_val_data[idx]);
        }
        return stat.calc(fill);
    }

    template <typename TOut>
    void operator()(const T* cond_row, const T* val_row, TOut* output) {
        ++m_count;
        for (int ii = 0; ii < m_column_size; ++ii) {
            auto* _start_cond_data = m_cond_data_.data() + ii * window_size;
            auto* _start_val_data = m_val_data_.data() + ii * window_size;

            output[ii] = handle(cond_row[ii], val_row[ii], _start_cond_data, _start_val_data);
        }
    }

    void set_param(const std::string& key, const std::string& value) {
        if (key == "method") {
            method = std::stoi(value);
        } else if (key == "quantile") {
            q = std::stod(value);
        } else if (key == "least") {
            least = std::stoi(value);
        } else if (key == "fill") {
            fill = std::stod(value);
        }
    }
};

using TsLteMean = rolling_cond_stat_rb_range<double, std::less_equal<double>, CondMean>;
using TsGteMean = rolling_cond_stat_rb_range<double, std::greater_equal<double>, CondMean>;
using TsLteMax = rolling_cond_stat_rb_range<double, std::less_equal<double>, CondMax>;
using TsGteMax = rolling_cond_stat_rb_range<double, std::greater_equal<double>, CondMax>;
using TsLteMin = rolling_cond_stat_rb_range<double, std::less_equal<double>, CondMin>;
using TsGteMin = rolling_cond_stat_rb_range<double, std::greater_equal<double>, CondMin>;
using TsLteSd = rolling_cond_stat_rb_range<double, std::less_equal<double>, CondSd>;
using TsGteSd = rolling_cond_stat_rb_range<double, std::greater_equal<double>, CondSd>;
}  // namespace ornate

#endif
