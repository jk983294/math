#ifndef ORNATE_MATH_CONTAINER_H
#define ORNATE_MATH_CONTAINER_H

#include "math_utils.h"
#include <math_vector.h>

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
        m_container.resize(window_size, std::vector<T>(m_column_size, T()));
    }

    void clear() {
        m_head_index = -1;
        m_count = 0;
    }

    void set(int size, int m_column_size_) {
        window_size = size + 1;
        m_column_size = m_column_size_;
        m_container.resize(window_size, std::vector<T>(m_column_size, T()));
    }

    const T *get_old_row() {
        return get_old_row_(m_count, m_head_index);
    }
    const T* get_next_old_row() {
        int next_idx = m_head_index + 1;
        if (next_idx == window_size) next_idx = 0;
        return get_old_row_(m_count + 1, next_idx);
    }

    const T *get_new_row() { return m_container[m_head_index].data(); }
    const T* get_next_new_row() {
        int next_idx = m_head_index + 1;
        if (next_idx == window_size) next_idx = 0;
        return m_container[next_idx].data();
    }

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

    double warm_up(bool all=false) {
        double ret = 0.;
        if (!all) {
            auto* ptr = get_next_old_row();
            if (ptr) ret += mean(ptr, m_column_size);
            auto* ptr1 = get_next_new_row();
            if (ptr1) ret += mean(ptr1, m_column_size);
        } else {
            for (const auto& d : m_container) {
                ret += mean(d);
            }
        }
        return ret;
    }

private:
    const T* get_old_row_(int count, int next_idx) {
        if (count >= window_size) {
            int old_index = next_idx - window_size + 1;
            if (old_index < 0) old_index += window_size;
            return m_container[old_index].data();
        } else {
            return nullptr;
        }
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

    void clear() {
        m_head_index = -1;
        m_count = 0;
    }

    void set(int size, int m_column_size_) {
        window_size = size + 1;
        m_column_size = m_column_size_;
        m_container.resize(window_size, nullptr);
    }

    const T *get_old_row() {
        return get_old_row_(m_count, m_head_index);
    }

    const T* get_next_old_row() {
        int next_idx = m_head_index + 1;
        if (next_idx == window_size) next_idx = 0;
        return get_old_row_(m_count + 1, m_head_index);
    }

    const T *get_ld_old_row() {
        if (m_count >= window_size) {
            int old_index = m_head_index - window_size + 1;
            if (old_index < 0) old_index += window_size;
            return m_container[old_index];
        } else {
            return m_container[0];
        }
    }

    const T *get_new_row() { return m_container[m_head_index]; }
    const T* get_next_new_row() {
        int next_idx = m_head_index + 1;
        if (next_idx == window_size) next_idx = 0;
        return m_container[next_idx];
    }

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

    double warm_up(bool all=false) {
        double ret = 0.;
        if (!all) {
            auto* ptr = get_next_old_row();
            if (ptr) ret += mean(ptr, m_column_size);
            auto* ptr1 = get_next_new_row();
            if (ptr1) ret += mean(ptr1, m_column_size);
        } else {
            for (const auto& d : m_container) {
                if (d) ret += mean(d, m_column_size);
            }
        }
        return ret;
    }

private:
    const T* get_old_row_(int count, int next_idx) {
        if (count >= window_size) {
            int old_index = next_idx - window_size + 1;
            if (old_index < 0) old_index += window_size;
            return m_container[old_index];
        } else {
            return nullptr;
        }
    }
};

template <typename T = double>
struct rolling_data_container2 {
    int m_size1{-1}, m_size2{-1};
    int window_size{-1};
    int m_column_size{-1};
    int m_head_index{-1};
    int m_count{0};
    std::vector<std::vector<T>> m_container;

    rolling_data_container2() = default;

    rolling_data_container2(int size1, int size2, int column_size_) { set(size1, size2, column_size_); }

    void clear() {
        m_head_index = -1;
        m_count = 0;
    }

    void set(int size1, int size2, int m_column_size_) {
        m_size1 = size1;
        m_size2 = size2;
        window_size = std::max(size1, size2) + 1;
        m_column_size = m_column_size_;
        m_container.resize(window_size);
        for (auto& c : m_container) c.resize(m_column_size, T());
    }

    const T *get_old_row0() {
        return get_old_row_(m_count, m_head_index, m_size1);
    }

    const T *get_old_row1() {
        return get_old_row_(m_count, m_head_index, m_size2);
    }

    const T* get_next_old_row0() {
        int next_idx = m_head_index + 1;
        if (next_idx == window_size) next_idx = 0;
        return get_old_row_(m_count + 1, next_idx, m_size1);
    }
    const T* get_next_old_row1() {
        int next_idx = m_head_index + 1;
        if (next_idx == window_size) next_idx = 0;
        return get_old_row_(m_count + 1, next_idx, m_size2);
    }

    const T *get_new_row() { return m_container[m_head_index].data(); }
    const T* get_next_new_row() {
        int next_idx = m_head_index + 1;
        if (next_idx == window_size) next_idx = 0;
        return m_container[next_idx].data();
    }

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

    double warm_up() {
        double ret = 0.;
        auto* ptr = get_next_old_row0();
        if (ptr) ret += mean(ptr, m_column_size);
        auto* ptr1 = get_next_new_row();
        if (ptr1) ret += mean(ptr1, m_column_size);
        auto* ptr2 = get_next_old_row1();
        if (ptr2) ret += mean(ptr2, m_column_size);
        return ret;
    }

private:
    const T* get_old_row_(int count, int next_idx, int size) {
        if (count >= size + 1) {
            int old_index = next_idx - size;
            if (old_index < 0) old_index += window_size;
            return m_container[old_index].data();
        } else {
            return nullptr;
        }
    }
};

template <typename T>
struct rolling_pointer_container2 {
    int m_size1{-1}, m_size2{-1};
    int window_size{-1};
    int m_column_size{-1};
    int m_head_index{-1};
    int m_count{0};
    std::vector<const T *> m_container;

    rolling_pointer_container2() = default;

    rolling_pointer_container2(int size1, int size2, int column_size_) { set(size1, size2, column_size_); }

    void clear() {
        m_head_index = -1;
        m_count = 0;
    }

    void set(int size1, int size2, int m_column_size_) {
        m_size1 = size1;
        m_size2 = size2;
        window_size = std::max(size1, size2) + 1;
        m_column_size = m_column_size_;
        m_container.resize(window_size, nullptr);
    }

    const T *get_old_row0() {
        return get_old_row_(m_count, m_head_index, m_size1);
    }

    const T *get_old_row1() {
        return get_old_row_(m_count, m_head_index, m_size2);
    }

    const T* get_next_old_row0() {
        int next_idx = m_head_index + 1;
        if (next_idx == window_size) next_idx = 0;
        return get_old_row_(m_count + 1, next_idx, m_size1);
    }
    const T* get_next_old_row1() {
        int next_idx = m_head_index + 1;
        if (next_idx == window_size) next_idx = 0;
        return get_old_row_(m_count + 1, next_idx, m_size2);
    }

    const T *get_new_row() { return m_container[m_head_index]; }
    const T* get_next_new_row() {
        int next_idx = m_head_index + 1;
        if (next_idx == window_size) next_idx = 0;
        return m_container[next_idx];
    }

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

    double warm_up() {
        double ret = 0.;
        auto* ptr = get_next_old_row0();
        if (ptr) ret += mean(ptr, m_column_size);
        auto* ptr1 = get_next_new_row();
        if (ptr1) ret += mean(ptr1, m_column_size);
        auto* ptr2 = get_next_old_row1();
        if (ptr2) ret += mean(ptr2, m_column_size);
        return ret;
    }

private:
    const T* get_old_row_(int count, int next_idx, int size) {
        if (count >= size + 1) {
            int old_index = next_idx - size;
            if (old_index < 0) old_index += window_size;
            return m_container[old_index];
        } else {
            return nullptr;
        }
    }
};
}  // namespace ornate

#endif
