#ifndef ORNATE_MATH_CONTAINER_H
#define ORNATE_MATH_CONTAINER_H

#include "math_utils.h"

namespace ornate {

namespace detail {
template <typename T, typename U>
void _data_copy2vector(const U *pData, T *_container, int size) {
    for (int ii = 0; ii < size; ++ii) {
        _container[ii] = pData[ii];
    }
}

template <typename T>
void _data_copy2vector(const T *pData, T *_container, int size) {
    std::copy(pData, pData + size, _container);
}

template <typename T, typename U>
void _data_copy2vector(const U *pData, std::vector<T> &_container, int size) {
    _data_copy2vector(pData, _container.data(), size);
}

template <typename T>
void _data_copy2vector(const T *pData, std::vector<T> &_container, int size) {
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
        m_container.resize(window_size);
        for (auto &c : m_container) c.resize(m_column_size);
    }

    const T *get_old_row() {
        if (m_count >= window_size) {
            int old_index = m_head_index - window_size + 1;
            if (old_index < 0) old_index += window_size;
            return m_container[old_index].data();
        } else {
            return nullptr;
        }
    }

    const T *get_new_row() { return m_container[m_head_index].data(); }

    /**
     * @param idx [0, ts_window)
     */
    const T *get_row_by_idx(int idx) {
        idx = m_head_index - idx;
        if (idx < 0) idx += window_size;
        return m_container[idx].data();
    }

    template <typename U>
    void push(const std::vector<U> &data) {
        ++m_count;
        ++m_head_index;
        if (m_head_index == window_size) m_head_index = 0;
        detail::_data_copy2vector(data.data(), m_container[m_head_index], m_column_size);
    }

    std::vector<T> &get_push_row() {
        ++m_count;
        ++m_head_index;
        if (m_head_index == window_size) m_head_index = 0;
        return m_container[m_head_index];
    }
};

template <typename T>
struct rolling_pointer_container {
    int window_size{-1};
    int m_column_size{-1};
    int m_head_index{-1};
    int m_count{0};
    std::vector<const T *> m_container;

    rolling_pointer_container() = default;

    rolling_pointer_container(int size, int m_column_size_) : window_size{size + 1}, m_column_size{m_column_size_} {
        m_container.resize(window_size, nullptr);
    }

    void set(int size, int m_column_size_) {
        window_size = size + 1;
        m_column_size = m_column_size_;
        m_container.resize(window_size, nullptr);
    }

    const T *get_old_row() {
        if (m_count >= window_size) {
            int old_index = m_head_index - window_size + 1;
            if (old_index < 0) old_index += window_size;
            return m_container[old_index];
        } else {
            return nullptr;
        }
    }

    const T *get_new_row() { return m_container[m_head_index]; }

    /**
     * @param idx [0, ts_window)
     */
    const T *get_row_by_idx(int idx) {
        idx = m_head_index - idx;
        if (idx < 0) idx += window_size;
        return m_container[idx];
    }

    void push(const T *data) {
        ++m_count;
        ++m_head_index;
        if (m_head_index == window_size) m_head_index = 0;
        m_container[m_head_index] = data;
    }
};
}  // namespace ornate

#endif
