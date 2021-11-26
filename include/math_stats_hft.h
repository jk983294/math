#ifndef ORNATE_MATH_STATS_HFT_H
#define ORNATE_MATH_STATS_HFT_H

#include <unordered_map>
#include <unordered_set>
#include "math_random.h"
#include "math_stats.h"
#include "math_utils.h"
#include "math_vector.h"

using std::isfinite;

namespace ornate {
template <typename T = double>
std::pair<size_t, size_t> hft_calc_na_count(const std::vector<std::vector<T>>& data_) {
    size_t cnt = 0, total = 0;
    for (const auto& i : data_) {
        cnt += calc_na_count(i.data(), i.size());
        total += i.size();
    }
    return {cnt, total};
}

template <typename T>
std::pair<double, double> hft_math_sign_ratio(const std::vector<std::vector<T>>& data_) {
    size_t neg_cnt_cum = 0, pos_cnt_cum = 0;
    for (const auto& i : data_) {
        int neg_cnt, pos_cnt;
        std::tie(neg_cnt, pos_cnt) = math_sign_count(i.data(), i.size());
        neg_cnt_cum += neg_cnt;
        pos_cnt_cum += pos_cnt;
    }

    if (neg_cnt_cum + pos_cnt_cum > 0) {
        double total_cnt = neg_cnt_cum + pos_cnt_cum;
        return {neg_cnt_cum / total_cnt, pos_cnt_cum / total_cnt};
    }
    return {NAN, NAN};
}

}  // namespace ornate

#endif
