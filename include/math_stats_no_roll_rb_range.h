#ifndef ORNATE_MATH_STATS_NO_ROLL_RB_RANGE_H
#define ORNATE_MATH_STATS_NO_ROLL_RB_RANGE_H

#include <algorithm>
#include <functional>
#include "math_container.h"

namespace ornate {
struct no_roll_mean_rb_range {
    struct stat {
        int m_valid_count{0};
        void clear() { m_valid_count = 0; }
    };
    int m_column_size;
    std::vector<stat> stats;

    explicit no_roll_mean_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    template <typename TOut>
    void init(TOut* output) {
        for (auto& stat : stats) stat.clear();
        std::fill(output, output + m_column_size, 0);
    }

    template <typename T, typename TOut>
    void operator()(int idx, const T* _row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            if (std::isfinite(_row[i])) {
                ++stats[i].m_valid_count;
                output[i] += _row[i];
            }
        }
    }

    template <typename TOut>
    void final_result(TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            if (stats[i].m_valid_count > 0)
                output[i] = output[i] / stats[i].m_valid_count;
            else
                output[i] = NAN;
        }
    }

    void set_row_size(int row) {}
};

struct no_roll_sum_rb_range {
    int m_column_size;

    explicit no_roll_sum_rb_range(int column_size_) : m_column_size{column_size_} {}

    template <typename TOut>
    void init(TOut* output) {
        std::fill(output, output + m_column_size, NAN);
    }

    template <typename T, typename TOut>
    void operator()(int idx, const T* _row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            if (!std::isfinite(output[i])) {
                output[i] = _row[i];
            } else if (std::isfinite(_row[i])) {
                output[i] += _row[i];
            }
        }
    }

    template <typename TOut>
    void final_result(TOut* output) {}

    void set_row_size(int row) {}
};

struct no_roll_decay_rb_range {
    struct stat {
        int m_valid_count{0};
        void clear() { m_valid_count = 0; }
    };
    int m_column_size;
    int m_row_size{0};  // window size
    std::vector<stat> stats;

    explicit no_roll_decay_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    void set_row_size(int row) { m_row_size = row; }

    template <typename TOut>
    void init(TOut* output) {
        for (auto& stat : stats) stat.clear();
        std::fill(output, output + m_column_size, 0);
    }

    template <typename T, typename TOut>
    void operator()(int idx, const T* _row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            if (std::isfinite(_row[i])) {
                stats[i].m_valid_count += (m_row_size - idx);
                output[i] += _row[i] * (m_row_size - idx);
            }
        }
    }

    template <typename TOut>
    void final_result(TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            if (stats[i].m_valid_count > 0)
                output[i] = output[i] / stats[i].m_valid_count;
            else
                output[i] = NAN;
        }
    }
};

struct no_roll_variance_rb_range {
    struct stat {
        double total_square_sum{0};
        int m_valid_count{0};
        void clear() {
            total_square_sum = 0;
            m_valid_count = 0;
        }
    };
    int m_column_size;
    std::vector<stat> stats;

    explicit no_roll_variance_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    template <typename TOut>
    void init(TOut* output) {
        for (auto& stat : stats) stat.clear();
        std::fill(output, output + m_column_size, 0);
    }

    template <typename T, typename TOut>
    void operator()(int idx, const T* _row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            if (std::isfinite(_row[i])) {
                ++stats[i].m_valid_count;
                output[i] += _row[i];  // total_sum
                stats[i].total_square_sum += _row[i] * _row[i];
            }
        }
    }

    template <typename TOut>
    void final_result(TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            if (stats[i].m_valid_count > 1) {
                double mean = output[i] / stats[i].m_valid_count;
                output[i] =
                    (stats[i].total_square_sum - mean * mean * stats[i].m_valid_count) / (stats[i].m_valid_count - 1);
            } else
                output[i] = NAN;
        }
    }

    void set_row_size(int row) {}
};

struct no_roll_std_rb_range {
    struct stat {
        double total_sum{0}, total_square_sum{0};
        int m_valid_count{0};
        void clear() {
            total_sum = total_square_sum = 0;
            m_valid_count = 0;
        }
    };
    int m_column_size;
    std::vector<stat> stats;

    explicit no_roll_std_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    template <typename TOut>
    void init(TOut* output) {
        for (auto& stat : stats) stat.clear();
        std::fill(output, output + m_column_size, 0);
    }

    template <typename T, typename TOut>
    void operator()(int idx, const T* _row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            if (std::isfinite(_row[i])) {
                ++stats[i].m_valid_count;
                stats[i].total_sum += _row[i];  // total_sum
                stats[i].total_square_sum += _row[i] * _row[i];
            }
        }
    }

    template <typename TOut>
    void final_result(TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            if (stats[i].m_valid_count > 1) {
                double mean = stats[i].total_sum / stats[i].m_valid_count;
                output[i] = std::sqrt((stats[i].total_square_sum - mean * mean * stats[i].m_valid_count) /
                                      (stats[i].m_valid_count - 1));
            } else
                output[i] = NAN;
        }
    }

    void set_row_size(int row) {}
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
    };
    int m_column_size;
    std::vector<stat> stats;

    explicit no_roll_zscore_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    template <typename TOut>
    void init(TOut* output) {
        for (auto& stat : stats) stat.clear();
    }

    template <typename T, typename TOut>
    void operator()(int idx, const T* _row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            if (std::isfinite(_row[i])) {
                ++stats[i].m_valid_count;
                stats[i].total_sum += _row[i];  // total_sum
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
            if (std::isfinite(output[i])) {
                double mean = stats[i].total_sum / stats[i].m_valid_count;
                double stddev = std::sqrt((stats[i].total_square_sum - mean * mean * stats[i].m_valid_count) /
                                          (stats[i].m_valid_count - 1));
                output[i] = (output[i] - mean) / stddev;
            }
        }
    }

    void set_row_size(int row) {}
};

struct no_roll_score_rb_range {
    struct stat {
        double total_sum{0};
        int m_valid_count{0};
        void clear() {
            total_sum = 0;
            m_valid_count = 0;
        }
    };
    int m_column_size;
    std::vector<stat> stats;

    explicit no_roll_score_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    template <typename TOut>
    void init(TOut* output) {
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
            if (std::isfinite(output[i])) {
                output[i] = output[i] - stats[i].total_sum / stats[i].m_valid_count;
            }
        }
    }

    void set_row_size(int row) {}
};

struct no_roll_cov_rb_range {
    struct stat {
        double sumy{0}, sumxy{0};
        int m_valid_count{0};
        void clear() {
            sumy = sumxy = 0;
            m_valid_count = 0;
        }
    };
    int m_column_size;
    std::vector<stat> stats;

    explicit no_roll_cov_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    template <typename TOut>
    void init(TOut* output) {
        for (auto& stat : stats) stat.clear();
        std::fill(output, output + m_column_size, 0);
    }

    template <typename T, typename TOut>
    void operator()(int idx, const T* x_row, const T* y_row, TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            if (std::isfinite(x_row[i]) && std::isfinite(y_row[i])) {
                ++stats[i].m_valid_count;
                stats[i].sumxy += x_row[i] * y_row[i];
                output[i] += x_row[i];  // sumx
                stats[i].sumy += y_row[i];
            }
        }
    }

    template <typename TOut>
    void final_result(TOut* output) {
        for (int i = 0; i < m_column_size; ++i) {
            if (stats[i].m_valid_count > 1) {
                double mean_x = output[i] / stats[i].m_valid_count;
                double mean_y = stats[i].sumy / stats[i].m_valid_count;
                output[i] = (stats[i].sumxy - mean_x * mean_y * stats[i].m_valid_count) / (stats[i].m_valid_count - 1);
            } else {
                output[i] = NAN;
            }
        }
    }

    void set_row_size(int row) {}
};

struct no_roll_corr_rb_range {
    struct stat {
        double sumx{0}, sumy{0}, sumxy{0}, sum_x2{0}, sum_y2{0};
        int m_valid_count{0};
        void clear() {
            sumx = sumy = sumxy = sum_x2 = sum_y2 = 0;
            m_valid_count = 0;
        }
    };
    int m_column_size;
    std::vector<stat> stats;

    explicit no_roll_corr_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    template <typename TOut>
    void init(TOut* output) {
        for (auto& stat : stats) stat.clear();
        std::fill(output, output + m_column_size, 0);
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
            if (stats[i].m_valid_count > 1) {
                double mean_x = stats[i].sumx / stats[i].m_valid_count;
                double mean_y = stats[i].sumy / stats[i].m_valid_count;
                double cov = (stats[i].sumxy - mean_x * mean_y * stats[i].m_valid_count) / (stats[i].m_valid_count - 1);
                double var1 =
                    (stats[i].sum_x2 - mean_x * mean_x * stats[i].m_valid_count) / (stats[i].m_valid_count - 1);
                double var2 =
                    (stats[i].sum_y2 - mean_y * mean_y * stats[i].m_valid_count) / (stats[i].m_valid_count - 1);
                if (std::isfinite(var1) && std::isfinite(var2)) {
                    double numerator = sqrt(var1) * sqrt(var2);
                    if (numerator < epsilon)
                        output[i] = NAN;
                    else
                        output[i] = cov / numerator;
                } else {
                    output[i] = NAN;
                }
            } else {
                output[i] = NAN;
            }
        }
    }

    void set_row_size(int row) {}
};

struct no_roll_skew_rb_range {
    struct stat {
        double total_x1{0}, total_x2{0}, total_x3{0};
        int m_valid_count{0};
        void clear() {
            total_x1 = total_x2 = total_x3 = 0;
            m_valid_count = 0;
        }
    };
    int m_column_size;
    std::vector<stat> stats;

    explicit no_roll_skew_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    template <typename TOut>
    void init(TOut* output) {
        for (auto& stat : stats) stat.clear();
        std::fill(output, output + m_column_size, NAN);
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
            if (stats[i].m_valid_count >= 2) {
                double mean = stats[i].total_x1 / stats[i].m_valid_count;
                double mean2 = mean * mean;
                double var = stats[i].total_x2 / stats[i].m_valid_count - mean2;
                if (var > 1e-14) {
                    double mean3 = mean * mean * mean;
                    double m3 = stats[i].total_x3 / stats[i].m_valid_count -
                                3 * mean * stats[i].total_x2 / stats[i].m_valid_count + 2 * mean3;
                    output[i] = m3 / std::pow(var, 1.5);
                }
            }
        }
    }

    void set_row_size(int row) {}
};

struct no_roll_kurtosis_rb_range {
    struct stat {
        double total_x1{0}, total_x2{0}, total_x3{0}, total_x4{0};
        ;
        int m_valid_count{0};
        void clear() {
            total_x1 = total_x2 = total_x3 = total_x4 = 0;
            m_valid_count = 0;
        }
    };
    int m_column_size;
    std::vector<stat> stats;

    explicit no_roll_kurtosis_rb_range(int column_size_) : m_column_size{column_size_} { stats.resize(m_column_size); }

    template <typename TOut>
    void init(TOut* output) {
        for (auto& stat : stats) stat.clear();
        std::fill(output, output + m_column_size, NAN);
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
            if (stats[i].m_valid_count >= 2) {
                double mean = stats[i].total_x1 / stats[i].m_valid_count;
                double mean2 = mean * mean;
                double var = stats[i].total_x2 / stats[i].m_valid_count - mean2;
                if (var > 1e-14) {
                    double mean3 = mean2 * mean;
                    double mean4 = mean3 * mean;
                    double m4 = stats[i].total_x4 / stats[i].m_valid_count -
                                4 * mean * stats[i].total_x3 / stats[i].m_valid_count +
                                6 * mean2 * stats[i].total_x2 / stats[i].m_valid_count - 3 * mean4;
                    output[i] = m4 / std::pow(var, 2) - 3.0;
                }
            }
        }
    }

    void set_row_size(int row) {}
};

}  // namespace ornate

#endif
