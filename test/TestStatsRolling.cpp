#include <math_utils.h>
#include <iostream>
#include "catch.hpp"
#include "math_stats.h"

using namespace std;
using namespace ornate;

static vector<double> x = {1, 2, NAN, 3, NAN, NAN, 4, 5.5, 6, 7, NAN, 8, 9};
static vector<double> data{1, 2, 3, 4, 5, 6};
static vector<double> data2{6, 5, 4, 3, 2, 1};
static vector<double> weight{1, 1, 1, 1, 1, 1};
static uint32_t small_window = 4;
static uint32_t shift_count = data.size() - small_window;

TEST_CASE("variance online", "[MathStatsRolling]") {
    double var = 0, mean = 0;
    size_t count = 1;
    for (auto obv : data) ornate::variance_online(var, mean, obv, count);
    REQUIRE(var == 3.5);
    REQUIRE(mean == 3.5);
    REQUIRE(count == 7);

    ornate::variance_rolling<> vr(data.size());
    for (auto obv : data) vr(obv);
    REQUIRE(vr.variance == 3.5);
    REQUIRE(vr.mean == 3.5);
    REQUIRE(vr.count == 6);
}

void test_variance_rolling(const vector<double>& _data, int window) {
    ornate::variance_rolling<> vr(window);
    rolling_variance_rb rvr(window);
    for (int i = 0; i < 6; ++i) {
        int to = i + 1;
        int from = to - window;
        if (from < 0) from = 0;
        double std_naive = ornate::std(_data.data() + from, to - from);
        double s = sqrtf(vr(_data[i]));
        rvr(_data[i]);
        if (std::isnan(std_naive)) {
            REQUIRE(std::isnan(s));
            REQUIRE(std::isnan(rvr.variance));
        } else {
            REQUIRE(sqrtf(vr.variance) == std_naive);
            REQUIRE(sqrtf(rvr.variance) == std_naive);
        }
    }
}

TEST_CASE("variance rolling", "[MathStatsRolling]") {
    vector<double> _data1{1, 2, 3, 4, 5, 6};
    test_variance_rolling(_data1, 3);
    test_variance_rolling(_data1, 5);

    vector<double> _data2{1, 2, NAN, 4, NAN, 6, 7, 8, 9, 10};
    test_variance_rolling(_data2, 3);
    test_variance_rolling(_data2, 5);

    vector<double> _data3{1.12, 1.13, NAN, 1.14, NAN, 1.12, 1.12, 1.13, 1.14, 1.12};
    test_variance_rolling(_data2, 3);
    test_variance_rolling(_data2, 5);
}

void test_covariance_rolling(const vector<double>& _data1, const vector<double>& _data2, int window) {
    ornate::covariance_rolling<> cvr(window);
    rolling_cov_rb rcr(window);
    for (int i = 0; i < 6; ++i) {
        int to = i + 1;
        int from = to - window;
        if (from < 0) from = 0;
        double cov_naive = ornate::cov(_data1.data() + from, _data2.data() + from, to - from);
        double s = cvr(_data1[i], _data2[i]);
        rcr(_data1[i], _data2[i]);
        if (std::isnan(cov_naive)) {
            REQUIRE(std::isnan(s));
            REQUIRE(std::isnan(rcr.covariance));
        } else {
            REQUIRE(FloatEqual(cvr.covariance, cov_naive));
            REQUIRE(FloatEqual(rcr.covariance, cov_naive));
        }
    }
}

TEST_CASE("covariance rolling", "[MathStatsRolling]") {
    test_covariance_rolling(data, data2, 3);
    test_covariance_rolling(data, data2, 5);

    vector<double> _data1{1, 2, NAN, 4, NAN, 6, 7, 8, 9, 10};
    vector<double> _data2{1.12, NAN, 1.13, 1.14, NAN, 1.12, 1.12, 1.13, 1.14, 1.12};
    test_covariance_rolling(_data1, _data2, 3);
    test_covariance_rolling(_data1, _data2, 5);
}

TEST_CASE("covariance online", "[MathStatsRolling]") {
    double covariance = 0, mean1 = 0, mean2 = 0;
    size_t count = 1;
    for (size_t i = 0; i < data.size(); ++i) {
        ornate::covariance_online(covariance, mean1, mean2, data[i], data2[i], count);
    }
    REQUIRE(covariance == -3.5);
    REQUIRE(mean1 == 3.5);
    REQUIRE(mean2 == 3.5);
    REQUIRE(count == 7);

    ornate::covariance_rolling<> cvr(data.size());
    for (size_t i = 0; i < data.size(); ++i) {
        cvr(data[i], data2[i]);
    }
    REQUIRE(cvr.covariance == -3.5);
    REQUIRE(cvr.meanA == 3.5);
    REQUIRE(cvr.meanB == 3.5);
    REQUIRE(cvr.count == 6);
}

TEST_CASE("covariance online small window", "[MathStatsRolling]") {
    double expected = -1.6666666666667;
    ornate::covariance_rolling<> cvr(small_window);
    for (size_t i = 0; i < data.size(); ++i) {
        cvr(data[i], data2[i]);
    }
    REQUIRE(FloatEqual(cvr.covariance, expected));
    REQUIRE(cvr.meanA == 4.5);
    REQUIRE(cvr.meanB == 2.5);
    REQUIRE(cvr.count == 4);

    REQUIRE(FloatEqual(covariance(data.data() + shift_count, data2.data() + shift_count, small_window), expected));
}

TEST_CASE("corr online", "[MathStatsRolling]") {
    REQUIRE(FloatEqual(ornate::corr(data, data2), -1.0));
    REQUIRE(FloatEqual(ornate::cov(data, data2), -3.5));

    OnlineCorrelation oc;
    for (size_t i = 0; i < data.size(); ++i) {
        oc.Push(data[i], data2[i]);
    }
    REQUIRE(FloatEqual(oc.Result(), -1.0));
}

void test_corr_rolling(const vector<double>& _data1, const vector<double>& _data2, int window) {
    ornate::corr_rolling cvr(window);
    rolling_corr_rb rcr(window);
    for (int i = 0; i < 6; ++i) {
        int to = i + 1;
        int from = to - window;
        if (from < 0) from = 0;
        double corr_naive = ornate::corr(_data1.data() + from, _data2.data() + from, to - from);
        double s = cvr(_data1[i], _data2[i]);
        rcr(_data1[i], _data2[i]);
        if (std::isnan(corr_naive)) {
            REQUIRE(std::isnan(s));
            REQUIRE(std::isnan(rcr.corr));
        } else {
            REQUIRE(FloatEqual(s, corr_naive));
            REQUIRE(FloatEqual(rcr.corr, corr_naive));
        }
    }
}

TEST_CASE("corr rolling", "[MathStatsRolling]") {
    corr_rolling cr(data.size());

    double ret = 0;
    for (size_t i = 0; i < data.size(); ++i) {
        ret = cr(data[i], data2[i]);
    }
    REQUIRE(FloatEqual(ret, -1.0));

    test_corr_rolling(data, data2, 3);
    test_corr_rolling(data, data2, 5);

    vector<double> _data1{1, 2, NAN, 4, NAN, 6, 7, 8, 9, 10};
    vector<double> _data2{1.12, NAN, 1.13, 1.14, NAN, 1.12, 1.12, 1.13, 1.14, 1.12};
    test_corr_rolling(_data1, _data2, 3);
    test_corr_rolling(_data1, _data2, 5);
}

static double calc_skew(const vector<double>& data_) {
    double mean_ = mean(data_);
    double std_ = 0;
    int valid_count = 0;
    for (double i : data_) {
        if (isfinite(i)) {
            ++valid_count;
            std_ += std::pow((i - mean_), 2);
        }
    }
    if (valid_count < 2) return NAN;
    std_ = sqrt(std_ / valid_count);
    if (std_ < 1e-7) return NAN;
    double ret = 0;
    for (double i : data_) {
        if (isfinite(i)) {
            ret += std::pow((i - mean_) / std_, 3);
        }
    }
    return ret / valid_count;
}

void test_skew_by_window(const vector<double>& x_, int window) {
    rolling_skew_rb sr(window);
    vector<double> y;

    double ret = 0;
    for (double d : x_) {
        if (y.size() < (size_t)window)
            y.push_back(d);
        else {
            for (int j = 1; j < window; ++j) {
                y[j - 1] = y[j];
            }
            y[window - 1] = d;
        }
        ret = sr(d);

        double expected = calc_skew(y);
        if (!FloatEqual(ret, expected)) {
            cerr << ret << " <-> " << expected << endl;
        }
        REQUIRE(FloatEqual(ret, expected));
    }
}

TEST_CASE("skew rolling nan", "[MathStatsRolling]") {
    test_skew_by_window(x, 3);
    test_skew_by_window(x, 4);
    test_skew_by_window(x, 5);
    test_skew_by_window(x, 6);
}

TEST_CASE("mean rolling", "[MathStatsRolling]") {
    mean_rolling<> sr(4);
    rolling_mean_rb rmr(4);

    double ret = 0;
    for (double i : data) {
        ret = sr(i);
        rmr(i);
        REQUIRE(FloatEqual(ret, rmr.mean));
    }
    REQUIRE(FloatEqual(ret, 4.5));
}

TEST_CASE("mean rolling nan", "[MathStatsRolling]") {
    mean_rolling<> sr(4);
    rolling_mean_rb rmr(4);

    double ret = 0;
    for (double i : x) {
        ret = sr(i);
        rmr(i);
        REQUIRE(FloatEqual(ret, rmr.mean));
    }
}

static double calc_kurtosis(const vector<double>& data_) {
    double mean_ = mean(data_);
    double std_ = 0;
    int valid_count = 0;
    for (double i : data_) {
        if (isfinite(i)) {
            ++valid_count;
            std_ += std::pow((i - mean_), 2);
        }
    }
    if (valid_count < 2) return NAN;
    std_ = sqrt(std_ / valid_count);
    if (std_ < 1e-7) return NAN;
    double ret = 0;
    for (double i : data_) {
        if (isfinite(i)) {
            ret += std::pow((i - mean_) / std_, 4);
        }
    }
    return ret / valid_count - 3.0;
}

void test_kurtosis_by_window(const vector<double>& x_, int window) {
    rolling_kurtosis_rb sr(window);
    vector<double> y;

    double ret = 0;
    for (double d : x_) {
        if (y.size() < (size_t)window)
            y.push_back(d);
        else {
            for (int j = 1; j < window; ++j) {
                y[j - 1] = y[j];
            }
            y[window - 1] = d;
        }
        ret = sr(d);

        double expected = calc_kurtosis(y);
        if (!FloatEqual(ret, expected)) {
            cerr << ret << " <-> " << expected << " window=" << window << " data=" << d << endl;
        }
        REQUIRE(FloatEqual(ret, expected));
    }
}

TEST_CASE("kurtosis rolling nan", "[MathStatsRolling]") {
    test_kurtosis_by_window(x, 3);
    test_kurtosis_by_window(x, 4);
    test_kurtosis_by_window(x, 5);
    test_kurtosis_by_window(x, 6);
}
