#ifndef ORNATE_MATH_DS_H
#define ORNATE_MATH_DS_H

#include <algorithm>
#include <cmath>
#include <vector>

namespace ornate {
struct DataSet {
    /**
     * align with add_row
     */
    void reserve(size_t n, size_t f_count) {
        m_n = n;
        m_f_count = f_count;
        m_y.reserve(n);
        m_xs.reserve(n * f_count);
    }

    bool add_row(double y, const std::vector<double>& xs) {
        if (!m_isRowMajor) return false;
        if (xs.size() != m_f_count) return false;
        m_y.push_back(y);
        std::copy(xs.begin(), xs.end(), std::back_inserter(m_xs));
        return true;
    }

    bool add_column_y(const std::vector<double>& y) {
        if (y.size() != m_n) return false;
        m_y = y;
        return true;
    }
    bool add_column_y_swap(std::vector<double>& y) {
        if (y.size() != m_n) return false;
        m_y.swap(y);
        return true;
    }

    bool add_column_x(const std::vector<double>& x) {
        if (m_isRowMajor) return false;
        if (x.size() != m_n) return false;
        std::copy(x.begin(), x.end(), std::back_inserter(m_xs));
        return true;
    }

    /**
     * align with add_column_x_to_row_major
     */
    void resize(size_t n, size_t f_count) {
        m_n = n;
        m_f_count = f_count;
        m_y.resize(n, NAN);
        m_xs.resize(n * f_count, NAN);
    }

    bool add_column_x_to_row_major(const std::vector<double>& x, int f_idx) {
        if (x.size() != m_n) return false;
        for (size_t i = 0; i < m_n; ++i) {
            m_xs[f_idx + i * m_f_count] = x[i];
        }
        return true;
    }

public:
    std::vector<double> m_y;
    std::vector<double> m_xs;
    size_t m_n{0};
    size_t m_f_count{0};
    bool m_isRowMajor{true};
};
}  // namespace ornate

#endif
