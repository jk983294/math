#include <math_utils.h>
#include <iostream>
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
    REQUIRE(ornate::mean_weighted(data.data(), weight.data(), data.size()) == 3.5);

    vector<float> x{3, 2, NAN, 1, NAN};
    REQUIRE(ornate::mean(x) == 2);
}

TEST_CASE("variance", "[MathStats]") {
    REQUIRE(variance(data.data(), data.size()) == 3.5);
    REQUIRE(variance_one_pass(data.data(), data.size()) == 3.5);
    REQUIRE(variance_two_pass(data.data(), data.size()) == 3.5);
    REQUIRE(variance_weighted(data.data(), weight.data(), data.size()) == 3.5);
}

TEST_CASE("variance online", "[MathStats]") {
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

TEST_CASE("covariance", "[MathStats]") {
    REQUIRE(ornate::covariance(data.data(), data2.data(), data.size()) == -3.5);
    REQUIRE(ornate::covariance_one_pass(data.data(), data2.data(), data.size()) == -3.5);
    REQUIRE(ornate::covariance_two_pass(data.data(), data2.data(), data.size()) == -3.5);
}

TEST_CASE("covariance online", "[MathStats]") {
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

TEST_CASE("covariance online small window", "[MathStats]") {
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

TEST_CASE("corr", "[MathStats]") {
    REQUIRE(FloatEqual(ornate::corr(data, data2), -1.0));
    REQUIRE(FloatEqual(ornate::cov(data, data2), -2.9166667));

    OnlineCorrelation oc;
    for (size_t i = 0; i < data.size(); ++i) {
        oc.Push(data[i], data2[i]);
    }
    REQUIRE(FloatEqual(oc.Result(), -1.0));
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
