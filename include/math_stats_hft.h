#ifndef ORNATE_MATH_STATS_HFT_H
#define ORNATE_MATH_STATS_HFT_H

#include <math_stats_once.h>
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
std::tuple<double, double, double> hft_math_sign_ratio(const std::vector<std::vector<T>>& data_) {
    rolling_sign_once rso;
    for (const auto& i : data_) {
        for (auto& d : i) rso(d);
    }
    return {rso.neg_ratio(), rso.full_neg_ratio(), rso.zero_ratio()};
}

template <typename T>
vector<double> hft_acf(const std::vector<std::vector<T>>& data_, int n_bar) {
    vector<double> ret(data_.size(), NAN);
    for (size_t ii = 0; ii < data_.size(); ++ii) {
        ret[ii] = ornate::acf(data_[ii], n_bar);
    }
    return ret;
}

template <typename T>
vector<double> hft_skew(const std::vector<std::vector<T>>& data_) {
    vector<double> ret(data_.size(), NAN);
    for (size_t ii = 0; ii < data_.size(); ++ii) {
        ret[ii] = ornate::math_skew(data_[ii]);
    }
    return ret;
}

template <typename T>
vector<double> hft_kurtosis(const std::vector<std::vector<T>>& data_) {
    vector<double> ret(data_.size(), NAN);
    for (size_t ii = 0; ii < data_.size(); ++ii) {
        ret[ii] = ornate::math_kurtosis(data_[ii]);
    }
    return ret;
}

template <typename T>
vector<double> hft_quantile(const std::vector<std::vector<T>>& data_, double q) {
    vector<double> ret(data_.size(), NAN);
    for (size_t ii = 0; ii < data_.size(); ++ii) {
        ret[ii] = ornate::quantile_const(data_[ii], q);
    }
    return ret;
}

template <typename T>
vector<double> hft_median(const std::vector<std::vector<T>>& data_) {
    vector<double> ret(data_.size(), NAN);
    for (size_t ii = 0; ii < data_.size(); ++ii) {
        ret[ii] = ornate::median(data_[ii]);
    }
    return ret;
}

template <typename T>
vector<double> hft_mean(const std::vector<std::vector<T>>& data_) {
    vector<double> ret(data_.size(), NAN);
    for (size_t ii = 0; ii < data_.size(); ++ii) {
        ret[ii] = ornate::mean(data_[ii]);
    }
    return ret;
}

template <typename T>
vector<double> hft_sd(const std::vector<std::vector<T>>& data_) {
    vector<double> ret(data_.size(), NAN);
    for (size_t ii = 0; ii < data_.size(); ++ii) {
        ret[ii] = ornate::std(data_[ii]);
    }
    return ret;
}

template <typename T>
vector<double> hft_max(const std::vector<std::vector<T>>& data_) {
    vector<double> ret(data_.size(), NAN);
    for (size_t ii = 0; ii < data_.size(); ++ii) {
        ret[ii] = ornate::v_max(data_[ii]);
    }
    return ret;
}

template <typename T>
vector<double> hft_min(const std::vector<std::vector<T>>& data_) {
    vector<double> ret(data_.size(), NAN);
    for (size_t ii = 0; ii < data_.size(); ++ii) {
        ret[ii] = ornate::v_min(data_[ii]);
    }
    return ret;
}

}  // namespace ornate

#endif
