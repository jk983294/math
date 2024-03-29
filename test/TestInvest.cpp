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
    REQUIRE(ornate::FloatEqual(nav1, nav[2]));
    REQUIRE(ornate::FloatEqual(nav1, nav[3]));
    REQUIRE(ornate::FloatEqual(nav4, nav[4]));
    REQUIRE(ornate::FloatEqual(nav5, nav[5]));
    REQUIRE(ornate::FloatEqual(nav6, nav[6]));

    nav = ornate::calc_bar_return_series(signals, rets, ins_num, true, 0, 0);
    REQUIRE(ornate::FloatEqual(nav0, nav[0]));
    REQUIRE(ornate::FloatEqual(nav1, nav[1]));
    REQUIRE(ornate::FloatEqual(nav1, nav[2]));
    REQUIRE(ornate::FloatEqual(nav1, nav[3]));
    REQUIRE(ornate::FloatEqual(nav4, nav[4]));
    REQUIRE(ornate::FloatEqual(nav5, nav[5]));
    REQUIRE(ornate::FloatEqual(nav6, nav[6]));
}

TEST_CASE("calc_bar_return_series 2 ins", "[calc_bar_return_series]") {
    std::vector<double> signals = {1, 2, NAN, NAN, -1, 2, -1, 2};
    std::vector<double> rets = {0.5, -0.5, NAN, NAN, -0.5, 0.5, 0, 0};
    int ins_num = 2;
    std::vector<double> nav = ornate::calc_bar_return_series(signals, rets, ins_num, false, 0, 0);
    double nav0 = ((1 + rets[0]) + (1 + rets[1])) / ins_num;
    double nav1 = nav0 * ((1 - rets[4]) + (1 + rets[5])) / ins_num;

    REQUIRE(ornate::FloatEqual(nav0, nav[0]));
    REQUIRE(ornate::FloatEqual(nav0, nav[1]));
    REQUIRE(ornate::FloatEqual(nav1, nav[2]));
    REQUIRE(ornate::FloatEqual(nav1, nav[3]));

    nav = ornate::calc_bar_return_series(signals, rets, ins_num, true, 0, 0);
    nav0 = 1.0 / 3. * (1 + rets[0]) + 2.0 / 3. * (1 + rets[1]);
    nav1 = nav0 * ((1 - rets[4]) + 2 * (1 + rets[5])) / 3.;
    REQUIRE(ornate::FloatEqual(nav0, nav[0]));
    REQUIRE(ornate::FloatEqual(nav0, nav[1]));
    REQUIRE(ornate::FloatEqual(nav1, nav[2]));
    REQUIRE(ornate::FloatEqual(nav1, nav[3]));
}

TEST_CASE("calc_bar_return_series 3 ins top 2", "[calc_bar_return_series]") {
    std::vector<double> signals = {1, 2, 0.5, NAN, NAN, NAN, -1, 2, 0.5, -1, 2, 0.5};
    std::vector<double> rets = {0.5, -0.5, -0.2, NAN, NAN, NAN, -0.5, 0.5, 0.2, 0, 0, 0};
    int ins_num = 3;
    int top_n = 2;
    std::vector<double> nav = ornate::calc_bar_return_series(signals, rets, ins_num, false, 0, 0, top_n);
    double nav0 = ((1 + rets[0]) + (1 + rets[1])) / top_n;
    double nav1 = nav0 * ((1 - rets[6]) + (1 + rets[7])) / top_n;

    REQUIRE(ornate::FloatEqual(nav0, nav[0]));
    REQUIRE(ornate::FloatEqual(nav0, nav[1]));
    REQUIRE(ornate::FloatEqual(nav1, nav[2]));
    REQUIRE(ornate::FloatEqual(nav1, nav[3]));

    nav = ornate::calc_bar_return_series(signals, rets, ins_num, true, 0, 0, top_n);
    nav0 = 1.0 / 3. * (1 + rets[0]) + 2.0 / 3. * (1 + rets[1]);
    nav1 = nav0 * ((1 - rets[6]) + 2 * (1 + rets[7])) / 3.;
    REQUIRE(ornate::FloatEqual(nav0, nav[0]));
    REQUIRE(ornate::FloatEqual(nav0, nav[1]));
    REQUIRE(ornate::FloatEqual(nav1, nav[2]));
    REQUIRE(ornate::FloatEqual(nav1, nav[3]));
}

TEST_CASE("calc_bar_return_series 3 ins top 2 sticky", "[calc_bar_return_series]") {
    std::vector<double> signals = {1, 2, 0.5, -1, 0.25, 0.5};
    std::vector<double> rets = {
        0.5, -0.5, -0.2, -0.5, 0.5, 0.2,
    };
    int ins_num = 3;
    int top_n = 2;
    std::vector<double> nav = ornate::calc_bar_return_series(signals, rets, ins_num, false, 0, 0, top_n, true);
    double nav0 = ((1 + rets[0]) + (1 + rets[1])) / top_n;
    double nav1 = nav0 * ((1 - rets[3]) + (1 + rets[4])) / 2;

    REQUIRE(ornate::FloatEqual(nav0, nav[0]));
    REQUIRE(ornate::FloatEqual(nav1, nav[1]));

    nav = ornate::calc_bar_return_series(signals, rets, ins_num, false, 0, 0, top_n, false);
    nav0 = ((1 + rets[0]) + (1 + rets[1])) / top_n;
    nav1 = nav0 * ((1 - rets[3]) + (1 + rets[5])) / 2.;
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
    REQUIRE(ornate::FloatEqual(nav1, nav[0][2]));
    REQUIRE(ornate::FloatEqual(nav1, nav[0][3]));
    REQUIRE(ornate::FloatEqual(nav4, nav[0][4]));
    REQUIRE(ornate::FloatEqual(nav5, nav[0][5]));
    REQUIRE(ornate::FloatEqual(nav6, nav[0][6]));

    REQUIRE(ornate::FloatEqual(nav0, nav[1][0]));
    REQUIRE(ornate::FloatEqual(nav1, nav[1][1]));
    REQUIRE(ornate::FloatEqual(nav1, nav[1][2]));
    REQUIRE(ornate::FloatEqual(nav1, nav[1][3]));
    REQUIRE(ornate::FloatEqual(nav4, nav[1][4]));
    REQUIRE(ornate::FloatEqual(nav5, nav[1][5]));
    REQUIRE(ornate::FloatEqual(nav6, nav[1][6]));
}

TEST_CASE("calc_max_dropdown", "[calc_max_dropdown]") {
    std::vector<double> signals = {1, 1.2, NAN, NAN, 0.9, 0.99, 1.2, 1.2};
    auto res = ornate::calc_max_dropdown_ratio(signals);
    REQUIRE(ornate::FloatEqual(res.first, 1 - (0.9 / 1.2)));
    REQUIRE(ornate::FloatEqual(res.second, double(4) / 8));
    signals.push_back(0.6);
    res = ornate::calc_max_dropdown_ratio(signals);
    REQUIRE(ornate::FloatEqual(res.first, 1 - (0.6 / 1.2)));
    REQUIRE(ornate::FloatEqual(res.second, double(8) / 9));
    signals.push_back(1.3);
    signals.push_back(0.6);
    res = ornate::calc_max_dropdown_ratio(signals);
    REQUIRE(ornate::FloatEqual(res.first, 1 - (0.6 / 1.3)));
    REQUIRE(ornate::FloatEqual(res.second, double(2) / 11));
}

TEST_CASE("clear last n ticks", "[calc_bar_return_series]") {
    std::vector<double> signals = {1, 1, -1, 1, -1, -1, 1, -1, 1, -1};
    std::vector<double> rets = {0.5, -0.5, -0.5, 0.5, -0.5, 0.5, -0.5, 0.5, -0.5, 0};
    TickSimStat stat(signals.size(), 1);
    stat.set_ti_num(signals.size());
    stat.set_clear_ticks({3, 4, 8, 9});
    std::vector<double> nav = stat.calc_bar_return_series(signals.data(), rets.data());
    double nav0 = 1 * (1 + rets[0]);
    double nav1 = nav0 * (1 + rets[1]);
    double nav2 = nav1 * (1 - rets[2]);
    double nav5 = nav2 * (1 - rets[5]);
    double nav6 = nav5 * (1 + rets[6]);
    double nav7 = nav6 * (1 - rets[7]);

    REQUIRE(ornate::FloatEqual(nav0, nav[0]));
    REQUIRE(ornate::FloatEqual(nav1, nav[1]));
    REQUIRE(ornate::FloatEqual(nav2, nav[2]));
    REQUIRE(ornate::FloatEqual(nav2, nav[3]));
    REQUIRE(ornate::FloatEqual(nav2, nav[4]));
    REQUIRE(ornate::FloatEqual(nav5, nav[5]));
    REQUIRE(ornate::FloatEqual(nav6, nav[6]));
    REQUIRE(ornate::FloatEqual(nav7, nav[7]));
    REQUIRE(ornate::FloatEqual(nav7, nav[8]));
    REQUIRE(ornate::FloatEqual(nav7, nav[9]));
}

TEST_CASE("calc_bar_return_series neutral 2 ins", "[calc_bar_return_series]") {
    std::vector<double> signals = {1, 2, NAN, NAN, -1, -2, -1, 2};
    std::vector<double> rets = {0.5, -0.5, NAN, NAN, -0.5, 0.5, 0, 0};
    int ins_num = 2, top_n = 2;
    std::vector<double> nav = ornate::calc_bar_return_series(signals, rets, ins_num, false, 0, 0, top_n, true, 0, true);
    double nav0 = ((1 - rets[0]) + (1 + rets[1])) / top_n;
    double nav1 = nav0 * ((1 + rets[4]) + (1 - rets[5])) / top_n;

    REQUIRE(ornate::FloatEqual(nav0, nav[0]));
    REQUIRE(ornate::FloatEqual(nav0, nav[1]));
    REQUIRE(ornate::FloatEqual(nav1, nav[2]));
    REQUIRE(ornate::FloatEqual(nav1, nav[3]));

    nav = ornate::calc_bar_return_series(signals, rets, ins_num, true, 0, 0, top_n, true, 0, true);
    REQUIRE(ornate::FloatEqual(nav0, nav[0]));
    REQUIRE(ornate::FloatEqual(nav0, nav[1]));
    REQUIRE(ornate::FloatEqual(nav1, nav[2]));
    REQUIRE(ornate::FloatEqual(nav1, nav[3]));
}

TEST_CASE("hft_calc_bar_return_series_vec 0", "[hft_calc_bar_return_series_vec]") {
    std::vector<double> signals = {1, 1, NAN, NAN, -1, 1, -1};
    std::vector<double> rets = {0.5, -0.5, NAN, NAN, -0.5, 0.5, 0};
    std::vector<double> nav = ornate::hft_calc_bar_return_series_vec(signals, rets, 0, 0);
    REQUIRE(nav.size() == 7);
    REQUIRE(nav == std::vector<double>({1, 1.5, 0.75, 0.75, 0.75, 0.75 * 1.5, 0.75 * 1.5 * 1.5}));
}

TEST_CASE("hft_calc_bar_return_series_vec 1", "[hft_calc_bar_return_series_vec]") {
    std::vector<double> signals = {1, 1, -1, 1, -1};
    std::vector<double> rets = {0.5, -0.5, -0.5, 0.5, 0};
    double cost = 0.1;
    std::vector<double> nav = ornate::hft_calc_bar_return_series_vec(signals, rets, 0, 0, cost);
    REQUIRE(nav.size() == 5);
    double nav1 = 1. * (1. + 0.5 - cost);
    double nav2 = 1. * ((1. + 0.5) * (1. - 0.5) - cost);
    double nav3 = nav2 * (1. + 0.5 - cost);
    double nav4 = nav3 * (1. + 0.5 - cost);
    REQUIRE(nav == std::vector<double>({1, nav1, nav2, nav3, nav4}));
}

TEST_CASE("hft_calc_bar_return_series_vec stop_ratio", "[hft_calc_bar_return_series_vec]") {
    std::vector<double> signals = {1, 1, 1, 1, -1};
    std::vector<double> rets = {-0.1, -0.1, -0.2, 0.5, 0};
    double cost = 0.;
    double stop_ratio = 0.3;
    std::vector<double> nav = ornate::hft_calc_bar_return_series_vec(signals, rets, 0, 0, cost, stop_ratio);
    double nav1 = 1. * (1. - 0.1 - cost);
    double nav2 = 1. * ((1. - 0.1) * (1. - 0.1) - cost);
    double nav3 = 1. * ((1. - 0.1) * (1. - 0.1) * (1. - 0.2) - cost);
    REQUIRE(nav == std::vector<double>({1, nav1, nav2, nav3, nav3}));
}

TEST_CASE("hft_calc_bar_return_series_vec stop_profit", "[hft_calc_bar_return_series_vec]") {
    std::vector<double> signals = {1, 1, 1, 1, -1, 1};
    std::vector<double> rets = {0.1, 0.1, 0.2, 0.5, 0.1, 0};
    double cost = 0.;
    double stop_profit = 0.3;
    std::vector<double> nav = ornate::hft_calc_bar_return_series_vec(signals, rets, 0, 0, cost, 0, stop_profit);
    double nav1 = 1. * (1. + 0.1 - cost);
    double nav2 = 1. * ((1. + 0.1) * (1. + 0.1) - cost);
    double nav3 = 1. * ((1. + 0.1) * (1. + 0.1) * (1. + 0.2) - cost);
    double nav4 = nav3 * (1. - 0.1 - cost);
    REQUIRE(nav == std::vector<double>({1, nav1, nav2, nav3, nav3, nav4}));
}

TEST_CASE("hft_calc_bar_return_series_vec stop_tick", "[hft_calc_bar_return_series_vec]") {
    std::vector<double> signals = {1, 1, 1, 1, -1, 1};
    std::vector<double> rets = {0.1, 0.1, 0.2, 0.5, 0.1, 0};
    double cost = 0.;
    int stop_tick = 3;
    std::vector<double> nav = ornate::hft_calc_bar_return_series_vec(signals, rets, 0, 0, cost, 0, 0, stop_tick);
    double nav1 = 1. * (1. + 0.1 - cost);
    double nav2 = 1. * ((1. + 0.1) * (1. + 0.1) - cost);
    double nav3 = 1. * ((1. + 0.1) * (1. + 0.1) * (1. + 0.2) - cost);
    double nav4 = nav3 * (1. - 0.1 - cost);
    REQUIRE(nav == std::vector<double>({1, nav1, nav2, nav3, nav3, nav4}));
}

TEST_CASE("TickRTStat stop_ratio", "[TickRTStat]") {
    std::vector<double> signals = {1, 1, 1, 1, -1};
    std::vector<double> rets = {-0.1, -0.1, -0.2, 0.5, 0};
    std::vector<double> close = {10, 10 * 0.9, 10 * 0.9 * 0.9, 10 * 0.9 * 0.9 * 0.8, 10 * 0.9 * 0.9 * 0.8 * 0.5};

    TickRTStat stat;
    stat.m_stop_ratio = 0.3;
    std::vector<int> pos;
    for (size_t i = 0; i < signals.size(); ++i) {
        pos.push_back(stat.process(signals[i], close[i]));
    }
    REQUIRE(pos == std::vector<int>({1, 1, 1, 0, -1}));
}

TEST_CASE("TickRTStat stop_profit", "[TickRTStat]") {
    std::vector<double> signals = {1, 1, 1, 1, -1, 1};
    std::vector<double> rets = {0.1, 0.1, 0.2, 0.5, 0.1, 0};
    std::vector<double> close = {
        10, 10 * 1.1, 10 * 1.1 * 1.1, 10 * 1.1 * 1.1 * 1.2, 10 * 1.1 * 1.1 * 1.2 * 0.5, 10 * 1.1 * 1.1 * 1.2 * 0.9};
    TickRTStat stat;
    stat.m_profit_ratio = 0.3;
    std::vector<int> pos;
    for (size_t i = 0; i < signals.size(); ++i) {
        pos.push_back(stat.process(signals[i], close[i]));
    }
    REQUIRE(pos == std::vector<int>({1, 1, 1, 0, -1, 1}));
}

TEST_CASE("TickRTStat stop_tick", "[TickRTStat]") {
    std::vector<double> signals = {1, 1, 1, 1, -1, 1};
    std::vector<double> rets = {0.1, 0.1, 0.2, 0.5, 0.1, 0};
    std::vector<double> close = {
        10, 10 * 1.1, 10 * 1.1 * 1.1, 10 * 1.1 * 1.1 * 1.2, 10 * 1.1 * 1.1 * 1.2 * 0.5, 10 * 1.1 * 1.1 * 1.2 * 0.9};

    TickRTStat stat;
    stat.m_stop_tick = 3;
    std::vector<int> pos;
    for (size_t i = 0; i < signals.size(); ++i) {
        pos.push_back(stat.process(signals[i], close[i]));
    }
    REQUIRE(pos == std::vector<int>({1, 1, 1, 0, -1, 1}));
}
