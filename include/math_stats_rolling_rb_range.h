#ifndef ORNATE_MATH_STATS_ROLLING_RB_RANGE_H
#define ORNATE_MATH_STATS_ROLLING_RB_RANGE_H

#include "math_utils.h"

namespace ornate {

template <typename T = double>
struct rolling_data_container {
    int window_size{-1};
    int m_column_size{-1};
    int m_head_index{-1};
    int m_count{0};
    std::vector<std::vector<T>> m_container;

    rolling_data_container() = default;
    rolling_data_container(int size, int m_column_size_) : window_size{size + 1}, m_column_size{m_column_size_} {
        m_container.resize(window_size, std::vector<T>(m_column_size));
    }

    void set(int size, int m_column_size_) {
        window_size = size + 1;
        m_column_size = m_column_size_;
        m_container.resize(window_size, std::vector<T>(m_column_size));
    }

    const T* get_old_row() {
        if (m_count >= window_size) {
            int old_index = m_head_index - window_size + 1;
            if (old_index < 0) old_index += window_size;
            return m_container[old_index].data();
        } else {
            return nullptr;
        }
    }

    const T* get_new_row() { return m_container[m_head_index].data(); }

    template <typename U>
    void push(const std::vector<U>& data) {
        ++m_count;
        ++m_head_index;
        if (m_head_index == window_size) m_head_index = 0;
        _data_copy2vector(data.data(), m_container[m_head_index]);
    }

private:
    template <typename U>
    void _data_copy2vector(const U* pData, std::vector<T>& _container) {
        for (int ii = 0; ii < m_column_size; ++ii) {
            _container[ii] = pData[ii];
        }
    }

    void _data_copy2vector(const T* pData, std::vector<T>& _container) {
        std::copy(pData, pData + m_column_size, _container.begin());
    }
};

template <typename T>
struct rolling_pointer_container {
    int window_size{-1};
    int m_column_size{-1};
    int m_head_index{-1};
    int m_count{0};
    std::vector<const T*> m_container;

    rolling_pointer_container() = default;
    rolling_pointer_container(int size, int m_column_size_) : window_size{size + 1}, m_column_size{m_column_size_} {
        m_container.resize(window_size, nullptr);
    }

    void set(int size, int m_column_size_) {
        window_size = size + 1;
        m_column_size = m_column_size_;
        m_container.resize(window_size, nullptr);
    }

    const T* get_old_row() {
        if (m_count >= window_size) {
            int old_index = m_head_index - window_size + 1;
            if (old_index < 0) old_index += window_size;
            return m_container[old_index];
        } else {
            return nullptr;
        }
    }

    const T* get_new_row() { return m_container[m_head_index]; }

    void push(const T* data) {
        ++m_count;
        ++m_head_index;
        if (m_head_index == window_size) m_head_index = 0;
        m_container[m_head_index] = data;
    }
};

struct rolling_mean_rb_range {
    struct stat {
        double total_sum{0};
        int m_valid_count{0};
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

            if (stat.m_valid_count > 0)
                output[i] = stat.total_sum / stat.m_valid_count;
            else
                output[i] = NAN;
        }
    }

    void set_row_size(int row) {}
};

struct rolling_sum_rb_range {
    struct stat {
        double total_sum{0};
        int m_valid_count{0};
    };
    int m_column_size;
    std::vector<stat> stats;

    explicit rolling_sum_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

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

            if (stat.m_valid_count > 0)
                output[i] = stat.total_sum;
            else
                output[i] = NAN;
        }
    }

    void set_row_size(int row) {}
};

struct rolling_variance_rb_range {
    struct stat {
        double total_sum{0}, total_square_sum{0};
        int m_valid_count{0};
    };
    int m_column_size;
    std::vector<stat> stats;

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

        if (st.m_valid_count > 1) {
            double mean = st.total_sum / st.m_valid_count;
            return (st.total_square_sum - mean * mean * st.m_valid_count) / (st.m_valid_count - 1);
        } else {
            return NAN;
        }
    }

    void set_row_size(int row) {}
};

struct rolling_cov_rb_range {
    struct stat {
        double sumx{0}, sumy{0}, sumxy{0};
        int m_valid_count{0};
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

        if (st.m_valid_count > 1) {
            double mean_x = st.sumx / st.m_valid_count;
            double mean_y = st.sumy / st.m_valid_count;
            return (st.sumxy - mean_x * mean_y * st.m_valid_count) / (st.m_valid_count - 1);
        } else {
            return NAN;
        }
    }

    void set_row_size(int row) {}
};

struct rolling_corr_rb_range {
    struct stat {
        double sumx{0}, sumy{0}, sumxy{0}, sum_x2{0}, sum_y2{0};
        int m_valid_count{0};
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

        if (st.m_valid_count > 1) {
            double mean_x = st.sumx / st.m_valid_count;
            double mean_y = st.sumy / st.m_valid_count;
            double cov = (st.sumxy - mean_x * mean_y * st.m_valid_count) / (st.m_valid_count - 1);
            double var1 = (st.sum_x2 - mean_x * mean_x * st.m_valid_count) / (st.m_valid_count - 1);
            double var2 = (st.sum_y2 - mean_y * mean_y * st.m_valid_count) / (st.m_valid_count - 1);
            if (std::isfinite(var1) && std::isfinite(var2)) {
                double numerator = sqrt(var1) * sqrt(var2);
                if (numerator < epsilon)
                    return NAN;
                else
                    return cov / numerator;
            }
        }
        return NAN;
    }

    void set_row_size(int row) {}
};

struct rolling_skew_rb_range {
    struct stat {
        double total_x1{0}, total_x2{0}, total_x3{0};
        int m_valid_count{0};
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

        if (st.m_valid_count >= 2) {
            double mean = st.total_x1 / st.m_valid_count;
            double mean2 = mean * mean;
            double var = st.total_x2 / st.m_valid_count - mean2;
            if (var <= 1e-14)
                return NAN;
            else {
                double mean3 = mean * mean * mean;
                double m3 = st.total_x3 / st.m_valid_count - 3 * mean * st.total_x2 / st.m_valid_count + 2 * mean3;
                return m3 / std::pow(var, 1.5);
            }
        }
        return NAN;
    }

    void set_row_size(int row) {}
};

struct rolling_kurtosis_rb_range {
    struct stat {
        double total_x1{0}, total_x2{0}, total_x3{0}, total_x4{0};
        int m_valid_count{0};
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

        if (st.m_valid_count >= 2) {
            double mean = st.total_x1 / st.m_valid_count;
            double mean2 = mean * mean;
            double var = st.total_x2 / st.m_valid_count - mean2;
            if (var <= 1e-14)
                return NAN;
            else {
                double mean3 = mean2 * mean;
                double mean4 = mean3 * mean;
                double m4 = st.total_x4 / st.m_valid_count - 4 * mean * st.total_x3 / st.m_valid_count +
                            6 * mean2 * st.total_x2 / st.m_valid_count - 3 * mean4;
                return m4 / std::pow(var, 2) - 3.0;
            }
        }
        return NAN;
    }

    void set_row_size(int row) {}
};

struct rolling_decay_rb_range {
    struct stat {
        double total_x1{0}, total{0};
        int m_valid_count{0}, m_valid_x1_count{0};
        ;
    };
    int m_column_size;
    int m_row_size{0};  // window size
    int count{0};
    std::vector<stat> stats;

    explicit rolling_decay_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    void set_row_size(int row) { m_row_size = row; }

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

        if (st.m_valid_count > 0) {
            return st.total / st.m_valid_count;
        }
        return NAN;
    }
};

}  // namespace ornate

#endif
