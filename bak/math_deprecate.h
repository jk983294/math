#ifndef ORNATE_MATH_DEPRECATE_H
#define ORNATE_MATH_DEPRECATE_H

#include <math_stats_rolling_rb.h>
#include <math_stats_rolling_rb_range.h>
#include <algorithm>
#include <functional>
#include "math_utils.h"

namespace ornate {
template <typename T>
struct rolling_rank_rb : public rolling_rank_base_rb<T> {
    explicit rolling_rank_rb(int size) : rolling_rank_base_rb<T>(size) {}

    using rolling_rank_base_rb<T>::m_valid_count;
    using rolling_rank_base_rb<T>::handle;

    double operator()(T new_value) {
        int lower, idx;
        std::tie(lower, idx) = handle(new_value);
        if (idx >= 0 && m_valid_count > 1)
            return (lower + idx) / (2.0 * (m_valid_count - 1));
        else
            return NAN;
    }
};

template <typename T>
struct rolling_rank_rb_range : public rolling_rank_base_rb_range<T> {
    using rolling_rank_base_rb_range<T>::m_count;
    using rolling_rank_base_rb_range<T>::m_column_size;
    using rolling_rank_base_rb_range<T>::stats;
    using rolling_rank_base_rb_range<T>::m_sorted_data_;
    using rolling_rank_base_rb_range<T>::window_size;
    using rolling_rank_base_rb_range<T>::handle;

    explicit rolling_rank_rb_range(int size) : rolling_rank_base_rb_range<T>(size) {}

    template <typename TOut>
    void operator()(const T* old_row, const T* new_row, TOut* output) {
        ++m_count;
        for (int i = 0; i < m_column_size; ++i) {
            auto& st = stats[i];
            auto* start_sorted_data = m_sorted_data_.data() + i * (window_size - 1);

            T old_value = get_nan<T>();
            if (old_row) old_value = old_row[i];
            int lower, idx, _valid_count;
            std::tie(lower, idx, _valid_count) = handle(new_row[i], old_value, st, start_sorted_data);
            if (idx >= 0 && _valid_count > 1)
                output[i] = (lower + idx) / (2.0 * (_valid_count - 1));
            else
                output[i] = NAN;
        }
    }

    void set_param(const std::string& key, const std::string& value) {}
};
}  // namespace ornate

#endif
