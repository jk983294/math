#ifndef ORNATE_MATH_STATS_NO_ROLL_RB_RANGE_H
#define ORNATE_MATH_STATS_NO_ROLL_RB_RANGE_H

#include <algorithm>
#include <functional>
#include "math_container.h"

namespace ornate {
struct no_roll_mean_rb_range {
    struct stat {
        double total_sum{0};
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

    explicit no_roll_mean_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    void init() {
        for (auto& stat : stats) stat.clear();
    }

    template <typename T, typename TOut>
    void operator()(int idx, const T* _row, TOut* output) {
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

struct no_roll_sum_rb_range {
    struct stat {
        double total_sum{0};
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

    explicit no_roll_sum_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    void init() {
        for (auto& stat : stats) stat.clear();
    }

    template <typename T, typename TOut>
    void operator()(int idx, const T* _row, TOut* output) {
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

struct no_roll_prod_rb_range {
    struct stat {
        double total_prod{1};
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

    explicit no_roll_prod_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    void init() {
        for (auto& stat : stats) stat.clear();
    }

    template <typename T, typename TOut>
    void operator()(int idx, const T* _row, TOut* output) {
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

struct no_roll_decay_rb_range {
    struct stat {
        double total{0};
        int m_valid_count{0};
        void clear() {
            total = 0;
            m_valid_count = 0;
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
    std::vector<stat> stats;

    explicit no_roll_decay_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    void set_row_size(int row) { m_row_size = row; }
    void set_param(const std::string& key, const std::string& value) {}

    void init() {
        for (auto& stat : stats) stat.clear();
    }

    template <typename T, typename TOut>
    void operator()(int idx, const T* _row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            if (std::isfinite(_row[i])) {
                stats[i].m_valid_count += (m_row_size - idx);
                stats[i].total += _row[i] * (m_row_size - idx);
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

struct no_roll_variance_rb_range {
    struct stat {
        double total_sum{0}, total_square_sum{0};
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
                    return (total_square_sum) / (m_valid_count - 1);
                }
            }
            return NAN;
        }
    };
    int m_column_size;
    std::vector<stat> stats;
    bool demean{true};  // for some series, its mean is 0, like return, no need to - mean

    explicit no_roll_variance_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    void init() {
        for (auto& stat : stats) stat.clear();
    }

    template <typename T, typename TOut>
    void operator()(int idx, const T* _row, TOut* output) {
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

struct no_roll_std_rb_range {
    struct stat {
        double total_sum{0}, total_square_sum{0};
        int m_valid_count{0};
        void clear() {
            total_sum = total_square_sum = 0;
            m_valid_count = 0;
        }
        double calc(bool demean_) const {
            if (m_valid_count > 1) {
                if (demean_) {
                    double mean = total_sum / m_valid_count;
                    return std::sqrt((total_square_sum - mean * mean * m_valid_count) / (m_valid_count - 1));
                } else {
                    return std::sqrt((total_square_sum) / (m_valid_count - 1));
                }
            }
            return NAN;
        }
    };
    int m_column_size;
    std::vector<stat> stats;
    bool demean{true};  // for some series, its mean is 0, like return, no need to - mean

    explicit no_roll_std_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    void init() {
        for (auto& stat : stats) stat.clear();
    }

    template <typename T, typename TOut>
    void operator()(int idx, const T* _row, TOut* output) {
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

struct no_roll_zscore_rb_range {
    struct stat {
        double total_sum{0}, total_square_sum{0};
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

    explicit no_roll_zscore_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    void init() {
        for (auto& stat : stats) stat.clear();
    }

    template <typename T, typename TOut>
    void operator()(int idx, const T* _row, TOut* output) {
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

struct no_roll_score_rb_range {
    struct stat {
        double total_sum{0};
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

    explicit no_roll_score_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    void init() {
        for (auto& stat : stats) stat.clear();
    }

    template <typename T, typename TOut>
    void operator()(int idx, const T* _row, TOut* output) {
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

struct no_roll_cov_rb_range {
    struct stat {
        double sumx{0}, sumy{0}, sumxy{0};
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

    explicit no_roll_cov_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    void init() {
        for (auto& stat : stats) stat.clear();
    }

    template <typename T, typename TOut>
    void operator()(int idx, const T* x_row, const T* y_row, TOut* output) {
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

struct no_roll_corr_rb_range {
    struct stat {
        double sumx{0}, sumy{0}, sumxy{0}, sum_x2{0}, sum_y2{0};
        int m_valid_count{0};
        void clear() {
            sumx = sumy = sumxy = sum_x2 = sum_y2 = 0;
            m_valid_count = 0;
        }
        double calc() const {
            if (m_valid_count > 1) {
                double mean_x = sumx / m_valid_count;
                double mean_y = sumy / m_valid_count;
                double cov = (sumxy - mean_x * mean_y * m_valid_count) / (m_valid_count - 1);
                double var1 = (sum_x2 - mean_x * mean_x * m_valid_count) / (m_valid_count - 1);
                double var2 = (sum_y2 - mean_y * mean_y * m_valid_count) / (m_valid_count - 1);
                if (std::isfinite(var1) && std::isfinite(var2)) {
                    double numerator = sqrt(var1) * sqrt(var2);
                    if (numerator >= epsilon) return cov / numerator;
                }
            }
            return NAN;
        }
    };
    int m_column_size;
    std::vector<stat> stats;

    explicit no_roll_corr_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    void init() {
        for (auto& stat : stats) stat.clear();
    }

    template <typename T, typename TOut>
    void operator()(int idx, const T* x_row, const T* y_row, TOut* output) {
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

struct no_roll_skew_rb_range {
    struct stat {
        double total_x1{0}, total_x2{0}, total_x3{0};
        int m_valid_count{0};
        void clear() {
            total_x1 = total_x2 = total_x3 = 0;
            m_valid_count = 0;
        }
        double calc() const {
            if (m_valid_count >= 2) {
                double mean = total_x1 / m_valid_count;
                double mean2 = mean * mean;
                double var = total_x2 / m_valid_count - mean2;
                if (var > 1e-14) {
                    double mean3 = mean * mean * mean;
                    double m3 = total_x3 / m_valid_count - 3 * mean * total_x2 / m_valid_count + 2 * mean3;
                    return m3 / std::pow(var, 1.5);
                }
            }
            return NAN;
        }
    };
    int m_column_size;
    std::vector<stat> stats;

    explicit no_roll_skew_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    void init() {
        for (auto& stat : stats) stat.clear();
    }

    template <typename T, typename TOut>
    void operator()(int idx, const T* _row, TOut* output) {
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

struct no_roll_kurtosis_rb_range {
    struct stat {
        double total_x1{0}, total_x2{0}, total_x3{0}, total_x4{0};
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

    explicit no_roll_kurtosis_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    void init() {
        for (auto& stat : stats) stat.clear();
    }

    template <typename T, typename TOut>
    void operator()(int idx, const T* _row, TOut* output) {
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

struct no_roll_regression2_rb_range {
    struct stat {
        double sum_x{0}, sum_y{0}, sum_x2{0}, sum_xy{0};
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

    explicit no_roll_regression2_rb_range(int column_size_) : m_column_size{column_size_} {
        stats.resize(m_column_size);
    }

    void init() {
        for (auto& stat : stats) stat.clear();
    }

    template <typename T>
    void operator()(int idx, const T* y_row, const T* x_row) {
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

    void add_new(stat& st, double y, double x) {
        if (std::isfinite(y) && std::isfinite(x)) {
            st.sum_x += x;
            st.sum_y += y;
            st.sum_x2 += x * x;
            st.sum_xy += x * y;
            ++st.m_valid_count;
        }
    }

    void set_row_size(int row) {}
    void set_param(const std::string& key, const std::string& value) {}
};

struct no_roll_regression3_rb_range {
    struct stat {
        double sum_x1_2{0}, sum_x2_2{0}, sum_y_2{0};
        double sum_x1{0}, sum_x2{0}, sum_y{0};
        double sum_x12{0}, sum_x1y{0}, sum_x2y{0};  // cross term
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

    explicit no_roll_regression3_rb_range(int column_size_) : m_column_size{column_size_} {
        stats.resize(m_column_size);
    }

    void init() {
        for (auto& stat : stats) stat.clear();
    }

    template <typename T>
    void operator()(int idx, const T* new_y, const T* new_x1, const T* new_x2) {
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

    void add_new(stat& st, double y, double x1, double x2) {
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

    void set_row_size(int row) {}
    void set_param(const std::string& key, const std::string& value) {}
};

struct no_roll_ema_rb_range {
    struct stat {
        double total_{0};
        int m_valid_count{0};
        void clear() {
            total_ = 0;
            m_valid_count = 0;
        }
        double calc() const {
            if (m_valid_count > 0)
                return total_;
            else
                return NAN;
        }
    };
    int m_column_size;
    double decay_coeff{0};
    std::vector<stat> stats;

    explicit no_roll_ema_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    void init() {
        for (auto& stat : stats) stat.clear();
    }

    template <typename T, typename TOut>
    void operator()(int idx, const T* _row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            if (std::isfinite(_row[i])) {
                ++stats[i].m_valid_count;
                if (stats[i].m_valid_count == 1)
                    stats[i].total_ = _row[i];
                else
                    stats[i].total_ = _row[i] * decay_coeff + stats[i].total_ * (1.0 - decay_coeff);
            } else {
                stats[i].total_ = stats[i].total_ * (1.0 - decay_coeff);
            }
        }
    }

    template <typename TOut>
    void final_result(TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            output[i] = stats[i].calc();
        }
    }

    void set_row_size(int ts_window) { decay_coeff = 2.0f / static_cast<double>(1 + ts_window); }
    void set_param(const std::string& key, const std::string& value) {}
};

struct no_roll_ema_hl_rb_range {
    struct stat {
        double total_{0}, total_w{0};
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
    double decay_coeff{0};
    std::vector<stat> stats;

    explicit no_roll_ema_hl_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    void init() {
        for (auto& stat : stats) stat.clear();
    }

    template <typename T, typename TOut>
    void operator()(int idx, const T* _row, TOut* output) {
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

    void set_row_size(int ts_window) { decay_coeff = ema_hl2decay(ts_window); }
    void set_param(const std::string& key, const std::string& value) {
        if (key == "hl" || key == "half_life") {
            decay_coeff = ema_hl2decay(std::stod(value));
        }
    }
};

}  // namespace ornate

#endif
