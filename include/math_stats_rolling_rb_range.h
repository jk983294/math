#ifndef ORNATE_MATH_STATS_ROLLING_RB_RANGE_H
#define ORNATE_MATH_STATS_ROLLING_RB_RANGE_H

#include <algorithm>
#include <functional>
#include "math_utils.h"

namespace ornate {

namespace detail {
template <typename T, typename U>
void _data_copy2vector(const U* pData, T* _container, int size) {
    for (int ii = 0; ii < size; ++ii) {
        _container[ii] = pData[ii];
    }
}

template <typename T>
void _data_copy2vector(const T* pData, T* _container, int size) {
    std::copy(pData, pData + size, _container);
}

template <typename T, typename U>
void _data_copy2vector(const U* pData, std::vector<T>& _container, int size) {
    _data_copy2vector(pData, _container.data(), size);
}

template <typename T>
void _data_copy2vector(const T* pData, std::vector<T>& _container, int size) {
    std::copy(pData, pData + size, _container.begin());
}
}  // namespace detail

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
        detail::_data_copy2vector(data.data(), m_container[m_head_index], m_column_size);
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
};

struct rolling_pct_rb_range {
    int m_column_size;

    explicit rolling_pct_rb_range(int column_size_) : m_column_size{column_size_} {}

    template <typename T, typename TOut>
    void operator()(const T* old_row, const T* new_row, TOut* output) {
        if (old_row) {
            for (int ii = 0; ii < m_column_size; ++ii) {
                output[ii] = new_row[ii] / old_row[ii] - 1.0f;
            }
        } else {
            std::fill(output, output + m_column_size, NAN);
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

struct rolling_std_rb_range {
    struct stat {
        double total_sum{0}, total_square_sum{0};
        int m_valid_count{0};
    };
    int m_column_size;
    std::vector<stat> stats;

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

        if (st.m_valid_count > 1) {
            double mean = st.total_sum / st.m_valid_count;
            return std::sqrt((st.total_square_sum - mean * mean * st.m_valid_count) / (st.m_valid_count - 1));
        } else {
            return NAN;
        }
    }

    void set_row_size(int row) {}
};

struct rolling_zscore_rb_range {
    struct stat {
        double total_sum{0}, total_square_sum{0};
        int m_valid_count{0};
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

    void set_row_size(int row) {}
};

struct rolling_score_rb_range {
    struct stat {
        double total_sum{0};
        int m_valid_count{0};
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

    float top_index(stat& st, Cell* start_cell) {
        if (st.front == st.rear) return NAN;
        return st.seq - start_cell[st.front].seq;
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
struct rolling_rank_rb_range {
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

    explicit rolling_rank_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }
    void set_row_size(int row) {
        window_size = row + 1;
        m_sorted_data_.resize(row * m_column_size);
    }

    template <typename TOut>
    void operator()(const T* old_row, const T* new_row, TOut* output) {
        ++m_count;
        for (int i = 0; i < m_column_size; ++i) {
            auto& st = stats[i];
            auto* start_sorted_data = m_sorted_data_.data() + i * (window_size - 1);

            T old_value = get_nan<T>();
            if (old_row) old_value = old_row[i];
            output[i] = calc(new_row[i], old_value, st, start_sorted_data);
        }
    }

    double calc(T new_value, T old_value, stat& st, SortItem* start_sorted_data) {
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

        if (idx >= 0 && st.m_valid_count > 1)
            return (lower + idx) / (2.0 * (st.m_valid_count - 1));
        else
            return NAN;
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
        while (pos < st.m_valid_count && start_sorted_data[pos].seq != seq && start_sorted_data[pos].data == seq) pos++;
        return pos;
    }
};

struct rolling_regression2_rb_range {
    struct stat {
        double sum_x{0}, sum_y{0}, sum_x2{0}, sum_xy{0};
        double a{NAN}, b{NAN};
        int m_valid_count{0};

        void calc() {
            if (m_valid_count > 1) {
                b = (m_valid_count * sum_xy - sum_x * sum_y) / (m_valid_count * sum_x2 - sum_x * sum_x);
                a = (sum_y - b * sum_x) / m_valid_count;
            } else {
                a = NAN;
                b = NAN;
            }
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
            a[i] = st.a;
            b[i] = st.b;
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

        st.calc();
    }

    void set_row_size(int row) {}
};

}  // namespace ornate

#endif
