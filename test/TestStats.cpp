#include <math_utils.h>
#include <iostream>
#include <numeric>
#include "catch.hpp"
#include "math_stats.h"

using namespace std;
using namespace ornate;

static vector<double> data{1, 2, 3, 4, 5, 6};
static vector<double> data2{6, 5, 4, 3, 2, 1};
static vector<double> weight{1, 1, 1, 1, 1, 1};
static uint32_t small_window = 4;
static uint32_t shift_count = data.size() - small_window;

TEST_CASE("mean", "[MathStats]") {
    REQUIRE(ornate::mean(data.data(), data.size()) == 3.5);
    REQUIRE(ornate::mean_weighted(data, weight) == 3.5);

    vector<float> x{3, 2, NAN, 1, NAN};
    REQUIRE(ornate::mean(x) == 2);
}

TEST_CASE("variance", "[MathStats]") {
    REQUIRE(variance(data.data(), data.size()) == 3.5);
    REQUIRE(variance_one_pass(data.data(), data.size()) == 3.5);
    REQUIRE(variance_two_pass(data.data(), data.size()) == 3.5);
    REQUIRE(variance_weighted(data.data(), weight.data(), data.size()) == 3.5);
}

TEST_CASE("covariance", "[MathStats]") {
    REQUIRE(ornate::covariance(data.data(), data2.data(), data.size()) == -3.5);
    REQUIRE(ornate::covariance_one_pass(data.data(), data2.data(), data.size()) == -3.5);
    REQUIRE(ornate::covariance_two_pass(data.data(), data2.data(), data.size()) == -3.5);
}

TEST_CASE("corr", "[MathStats]") {
    REQUIRE(FloatEqual(ornate::corr(data, data2), -1.0));
    REQUIRE(FloatEqual(ornate::cov(data, data2), -3.5));
}

TEST_CASE("std", "[MathStats]") {
    cout << setprecision(12) << ornate::std(data) << endl;
    REQUIRE(FloatEqual(ornate::std(data), 1.870828747));
}

TEST_CASE("regression", "[MathStats]") {
    double a, b, R;
    bool ret = ornate::regression(data, data2, &a, &b, &R);
    cout << setprecision(12) << a << " " << b << " " << R << endl;
    REQUIRE(ret);
    REQUIRE(FloatEqual(a, 7.0));
    REQUIRE(FloatEqual(b, -1.0));
    REQUIRE(FloatEqual(R, 1.0));
}

/**
 * y = -1 + 2 * x1 + 1 * x2
 */
TEST_CASE("regression3", "[MathStats]") {
    float b0, b1, b2;
    //    vector<float> x1 = {3, 4, 5, 6, 2};
    //    vector<float> x2 = {8, 5, 7, 3, 1};
    //    vector<float> y = {-3.7, 3.5, 2.5, 11.5, 5.7};

    vector<float> x1 = {0, 1, 2, 0, 0, 1, 1, 2, 2};
    vector<float> x2 = {0, 0, 0, 1, 2, 1., 2, 1, 2};
    vector<float> y(x1.size(), 0);
    vector<float> coeff = {-1, 2, 1};
    for (std::size_t i = 0; i < x1.size(); ++i) {
        y[i] = coeff[0] + coeff[1] * x1[i] + coeff[2] * x2[i];
    }
    bool ret = ornate::regression3(y, x1, x2, &b0, &b1, &b2);
    cout << setprecision(12) << b0 << " " << b1 << " " << b2 << endl;
    REQUIRE(ret);
    REQUIRE(FloatEqual(b0, coeff[0]));
    REQUIRE(FloatEqual(b1, coeff[1]));
    REQUIRE(FloatEqual(b2, coeff[2]));
}

TEST_CASE("ols", "[MathStats]") {
    float coeff = -2;
    vector<float> x1 = {0, 1, 2, 3, 4};
    vector<float> y(x1.size(), 0);
    for (std::size_t i = 0; i < x1.size(); ++i) {
        y[i] = coeff * x1[i];
    }
    float coeff_calc = ornate::ols(y, x1);
    REQUIRE(FloatEqual(coeff_calc, coeff));
}

/**
 * y = 2 * x1 + 1 * x2
 */
TEST_CASE("ols2", "[MathStats]") {
    float b1, b2;
    vector<float> x1 = {0, 1, 2, 0, 0, 1, 1, 2, 2};
    vector<float> x2 = {0, 0, 0, 1, 2, 1., 2, 1, 2};
    vector<float> y(x1.size(), 0);
    vector<float> coeff = {2, 1};
    for (std::size_t i = 0; i < x1.size(); ++i) {
        y[i] = coeff[0] * x1[i] + coeff[1] * x2[i];
    }
    bool ret = ornate::ols(y, x1, x2, &b1, &b2);
    cout << setprecision(12) << b1 << " " << b2 << endl;
    REQUIRE(ret);
    REQUIRE(FloatEqual(b1, coeff[0]));
    REQUIRE(FloatEqual(b2, coeff[1]));
}

TEST_CASE("quantile", "[MathStats]") {
    std::vector<double> v(10);
    std::iota(v.begin(), v.end(), 0.0);
    REQUIRE(FloatEqual(ornate::quantile(v, -0.1), 0));
    REQUIRE(ornate::quantile(v, 0.0) == 0);
    REQUIRE(FloatEqual(ornate::quantile(v, 0.1), 0.9));
    REQUIRE(FloatEqual(ornate::quantile(v, 0.8), 7.2));
    REQUIRE(FloatEqual(ornate::quantile(v, 0.9), 8.1));
    REQUIRE(ornate::quantile(v, 1.0) == 9);
    REQUIRE(ornate::quantile(v, 1.1) == 9);
}

TEST_CASE("ema_hl", "[MathStats]") {
    std::vector<double> v(10);
    std::iota(v.begin(), v.end(), 0.0);
    REQUIRE(FloatEqual(ornate::ema_hl(v, 5, 5, 4), 3.342143));
}

TEST_CASE("rank_last", "[MathStats]") {
    std::vector<double> v = {1, 2, 3, 4, 1};
    REQUIRE(FloatEqual(ornate::rank_last(v.data(), 4, 0), 0.125));
    REQUIRE(FloatEqual(ornate::rank_last(v.data(), 4, 5), 0.125));
    REQUIRE(FloatEqual(ornate::rank_last(v.data(), 3, 5), NAN));
    //    REQUIRE(FloatEqual(ornate::rank_last(v.data(), 4, 0), 0.25));
    //    REQUIRE(FloatEqual(ornate::rank_last(v.data(), 4, 5), 0.25));
}
