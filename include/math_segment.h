#ifndef ORNATE_MATH_SEGMENT_H
#define ORNATE_MATH_SEGMENT_H

#include <limits>
#include <vector>
#include "math_stats_rolling_rb_range.h"
#include "math_type.h"
#include "math_vector.h"

namespace ornate {
template <typename T>
std::vector<double> segment_sma(const T* signals, int total, int ins_num, int window) {
    std::vector<double> ret(total, NAN);
    rolling_mean_rb_range rmrr(ins_num);
    rolling_pointer_container<T> container(window, ins_num);
    for (int offset = 0; offset < total; offset += ins_num) {
        int valid_cnt = 0;
        for (int ii = 0; ii < ins_num; ++ii) {
            if (std::isfinite(signals[offset + ii])) ++valid_cnt;
        }

        if (valid_cnt > 0) {
            container.push(signals + offset);
            rmrr(container.get_old_row(), container.get_new_row(), ret.data() + offset);
        }
    }
    return ret;
}

template <typename T>
std::vector<T> segment_extract_skip_ti(const T* data, int total, int ins_num, int ti_num,
                                       const std::vector<int>& skip_ti) {
    std::vector<T> ret;
    ret.reserve(total);
    if (total % (ins_num * ti_num) == 0) {
        std::vector<bool> ti_index(ti_num, true);
        for (auto idx : skip_ti) ti_index[idx] = false;
        int di_num = total / (ins_num * ti_num);
        for (int di = 0; di < di_num; ++di) {
            for (int ti = 0; ti < ti_num; ++ti) {
                if (ti_index[ti]) {
                    int offset = di * ins_num * ti_num + ti * ins_num;
                    std::copy(data + offset, data + offset + ins_num, std::back_inserter(ret));
                }
            }
        }
    }
    return ret;
}

}  // namespace ornate

#endif
