#include <math_utils.h>
#include <iostream>
#include "catch.hpp"
#include "math_invest.h"

using namespace std;
using namespace ornate;

TEST_CASE("calc_bar_return_series 1 ins", "[calc_bar_return_series]") {
    std::vector<double> signals = {1, 1, NAN, NAN, -1, 1, -1};
    std::vector<double> rets = {0.5, -0.5, NAN, NAN, -0.5, 0.5, 0};
    int ins_num = 1;
    std::vector<double> nav = ornate::calc_bar_return_series(signals, rets, ins_num, false, 0, 0);
    double nav0 = 1 * (1 + rets[0]);
    double nav1 = nav0 * (1 + rets[1]);
    double nav4 = nav1 * (1 - rets[4]);
    double nav5 = nav4 * (1 + rets[5]);
    double nav6 = nav5 * (1 - rets[6]);
    REQUIRE(ornate::FloatEqual(nav0, 1.5));
    REQUIRE(ornate::FloatEqual(nav1, 1.5 * 0.5));
    REQUIRE(ornate::FloatEqual(nav4, 1.5 * 0.5 * 1.5));
    REQUIRE(ornate::FloatEqual(nav5, 1.5 * 0.5 * 1.5 * 1.5));
    REQUIRE(ornate::FloatEqual(nav6, 1.5 * 0.5 * 1.5 * 1.5 * 1.0));

    REQUIRE(ornate::FloatEqual(nav0, nav[0]));
    REQUIRE(ornate::FloatEqual(nav1, nav[1]));
    REQUIRE(ornate::FloatEqual(nav4, nav[2]));
    REQUIRE(ornate::FloatEqual(nav5, nav[3]));
    REQUIRE(ornate::FloatEqual(nav6, nav[4]));

    nav = ornate::calc_bar_return_series(signals, rets, ins_num, true, 0, 0);
    REQUIRE(ornate::FloatEqual(nav0, nav[0]));
    REQUIRE(ornate::FloatEqual(nav1, nav[1]));
    REQUIRE(ornate::FloatEqual(nav4, nav[2]));
    REQUIRE(ornate::FloatEqual(nav5, nav[3]));
    REQUIRE(ornate::FloatEqual(nav6, nav[4]));
}

TEST_CASE("calc_bar_return_series 2 ins", "[calc_bar_return_series]") {
    std::vector<double> signals = {1, 2, NAN, NAN, -1, 2, -1, 2};
    std::vector<double> rets = {0.5, -0.5, NAN, NAN, -0.5, 0.5, 0, 0};
    int ins_num = 2;
    std::vector<double> nav = ornate::calc_bar_return_series(signals, rets, ins_num, false, 0, 0);
    double nav0 = ((1 + rets[0]) + (1 + rets[1])) / ins_num;
    double nav1 = nav0 * ((1 - rets[4]) + (1 + rets[5])) / ins_num;

    REQUIRE(ornate::FloatEqual(nav0, nav[0]));
    REQUIRE(ornate::FloatEqual(nav1, nav[1]));

    nav = ornate::calc_bar_return_series(signals, rets, ins_num, true, 0, 0);
    nav0 = 1.0 / 3. * (1 + rets[0]) + 2.0 / 3. * (1 + rets[1]);
    nav1 = nav0 * ((1 - rets[4]) + 2 * (1 + rets[5])) / 3.;
    REQUIRE(ornate::FloatEqual(nav0, nav[0]));
    REQUIRE(ornate::FloatEqual(nav1, nav[1]));
}

TEST_CASE("calc_return_series_by_ii 2 ins", "[calc_return_series_by_ii]") {
    std::vector<double> signals = {1, 1, 1, 1, NAN, NAN, NAN, NAN, -1, -1, 1, 1, -1, -1};
    std::vector<double> rets = {0.5, 0.5, -0.5, -0.5, NAN, NAN, NAN, NAN, -0.5, -0.5, 0.5, 0.5, 0, 0};
    int ins_num = 2;
    std::vector<std::vector<double>> nav = ornate::calc_return_series_by_ii(signals, rets, ins_num, 0, 0);
    double nav0 = 1 * (1 + rets[0]);
    double nav1 = nav0 * (1 + rets[2]);
    double nav4 = nav1 * (1 - rets[8]);
    double nav5 = nav4 * (1 + rets[10]);
    double nav6 = nav5 * (1 - rets[12]);
    REQUIRE(ornate::FloatEqual(nav0, 1.5));
    REQUIRE(ornate::FloatEqual(nav1, 1.5 * 0.5));
    REQUIRE(ornate::FloatEqual(nav4, 1.5 * 0.5 * 1.5));
    REQUIRE(ornate::FloatEqual(nav5, 1.5 * 0.5 * 1.5 * 1.5));
    REQUIRE(ornate::FloatEqual(nav6, 1.5 * 0.5 * 1.5 * 1.5 * 1.0));

    REQUIRE(ornate::FloatEqual(nav0, nav[0][0]));
    REQUIRE(ornate::FloatEqual(nav1, nav[0][1]));
    REQUIRE(ornate::FloatEqual(nav4, nav[0][2]));
    REQUIRE(ornate::FloatEqual(nav5, nav[0][3]));
    REQUIRE(ornate::FloatEqual(nav6, nav[0][4]));

    REQUIRE(ornate::FloatEqual(nav0, nav[1][0]));
    REQUIRE(ornate::FloatEqual(nav1, nav[1][1]));
    REQUIRE(ornate::FloatEqual(nav4, nav[1][2]));
    REQUIRE(ornate::FloatEqual(nav5, nav[1][3]));
    REQUIRE(ornate::FloatEqual(nav6, nav[1][4]));
}

TEST_CASE("calc_max_dropdown", "[calc_max_dropdown]") {
    std::vector<double> signals = {1, 1.2, NAN, NAN, 0.9, 0.99, 1.2, 1.2};
    REQUIRE(ornate::FloatEqual(ornate::calc_max_dropdown_ratio(signals), 1 - (0.9 / 1.2)));
    signals.push_back(0.6);
    REQUIRE(ornate::FloatEqual(ornate::calc_max_dropdown_ratio(signals), 1 - (0.6 / 1.2)));
    signals.push_back(1.3);
    signals.push_back(0.65);
    REQUIRE(ornate::FloatEqual(ornate::calc_max_dropdown_ratio(signals), 1 - (0.65 / 1.3)));
}