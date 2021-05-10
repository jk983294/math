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

}  // namespace ornate

#endif
