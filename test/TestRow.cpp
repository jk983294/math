#include <math_row.h>
#include <algorithm>
#include <catch.hpp>

using namespace ornate;

TEST_CASE("row_sum", "[row_sum]") {
    std::vector<double> a = {1, 2, 3, 4, 5};
    auto res = row_sum(a.data(), a.size());
    REQUIRE(res == 15);
}

TEST_CASE("row_sum2", "[row_sum2]") {
    std::vector<double> a = {1, 2, 3, 4, 5};
    double expected = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        expected += a[i] * a[i];
    }
    auto res = row_sum2(a.data(), a.data(), a.size());
    REQUIRE(res == expected);
}

TEST_CASE("row_mean", "[row_mean]") {
    std::vector<double> a = {1, 2, 3, 4, 5};
    auto res = row_mean(a.data(), a.size());
    REQUIRE(res == 3);
}

TEST_CASE("row_mean2", "[row_mean2]") {
    std::vector<double> a = {1, 2, 3, 4, 5};
    double expected = 0, weight_ = 0;
    for (size_t i = 0; i < a.size(); ++i) {
        expected += a[i] * a[i];
        weight_ += a[i];
    }
    expected /= weight_;
    auto res = row_mean2(a.data(), a.data(), a.size());
    REQUIRE(res == expected);
}

TEST_CASE("row_cut", "[row_cut]") {
    std::vector<double> a = {1, 2, 3, 4, 5};
    auto res = row_cut(a.data(), 0.5, a.size());
    REQUIRE(res == 1);
    res = row_cut(a.data(), 1, a.size());
    REQUIRE(res == 1);
    res = row_cut(a.data(), 1.5, a.size());
    REQUIRE(res == 2);
    res = row_cut(a.data(), 2, a.size());
    REQUIRE(res == 2);
    res = row_cut(a.data(), 4.5, a.size());
    REQUIRE(res == 5);
    res = row_cut(a.data(), 5, a.size());
    REQUIRE(res == 5);
    res = row_cut(a.data(), 6, a.size());
    REQUIRE(res == 5);
}

TEST_CASE("row_sum_cut", "[row_sum_cut]") {
    std::vector<double> a = {1, 2, 3, 4, 5};
    auto res = row_sum_cut(a.data(), 0.5, a.size());
    REQUIRE(res == 1);
    res = row_sum_cut(a.data(), 1, a.size());
    REQUIRE(res == 1);
    res = row_sum_cut(a.data(), 1.5, a.size());
    REQUIRE(res == 2);
    res = row_sum_cut(a.data(), 2, a.size());
    REQUIRE(res == 2);
    res = row_sum_cut(a.data(), 14.5, a.size());
    REQUIRE(res == 5);
    res = row_sum_cut(a.data(), 15, a.size());
    REQUIRE(res == 5);
    res = row_sum_cut(a.data(), 16, a.size());
    REQUIRE(res == 5);
}

TEST_CASE("row_max", "[row_max]") {
    std::vector<double> a = {1, 2, 3, 4, 5};
    auto res = row_max(a.data(), a.size());
    REQUIRE(res == 5);
    auto idx = row_which_max(a.data(), a.size());
    REQUIRE(idx == 5);
}

TEST_CASE("row_min", "[row_min]") {
    std::vector<double> a = {1, 2, 3, 4, 5};
    auto res = row_min(a.data(), a.size());
    REQUIRE(res == 1);
    auto idx = row_which_min(a.data(), a.size());
    REQUIRE(idx == 1);
}

TEST_CASE("row_which_op", "[row_which_op]") {
    std::vector<double> a = {1, 2, 3, 4, 5};
    auto idx = row_which_eq(a.data(), 2, a.size());
    REQUIRE(idx == 2);
    REQUIRE(row_extract(a.data(), idx, a.size()) == 2);
    idx = row_which_gt(a.data(), 2, a.size());
    REQUIRE(idx == 3);
    REQUIRE(row_extract(a.data(), idx, a.size()) == 3);
    idx = row_which_gte(a.data(), 2, a.size());
    REQUIRE(idx == 2);
    REQUIRE(row_extract(a.data(), idx, a.size()) == 2);
    idx = row_which_lt(a.data(), 2, a.size());
    REQUIRE(idx == 1);
    REQUIRE(row_extract(a.data(), idx, a.size()) == 1);
    idx = row_which_lte(a.data(), 2, a.size());
    REQUIRE(idx == 1);
    REQUIRE(row_extract(a.data(), idx, a.size()) == 1);
}

TEST_CASE("row_sd", "[row_sd]") {
    std::vector<double> a = {1, 2, 3, 4, 5};
    auto res = row_sd(a.data(), a.size());
    REQUIRE(res == Approx(1.5811388301));
}

TEST_CASE("row_sd2", "[row_sd2]") {
    std::vector<double> a = {1, 2, 3, 4, 5};
    auto res = row_sd2(a.data(), a.data(), a.size());
    REQUIRE(res == Approx(1.3944333776));
}

TEST_CASE("row_median", "[row_median]") {
    std::vector<double> a = {1, 2, 3, 4, 5};
    auto res = row_median(a.data(), a.size());
    REQUIRE(res == Approx(3));
}

TEST_CASE("row_slope", "[row_slope]") {
    std::vector<double> a = {1, 2, 3, 4, 5};
    auto res = row_slope(a.data(), a.size());
    REQUIRE(res == Approx(1));
}

TEST_CASE("row_slope2", "[row_slope2]") {
    std::vector<double> a = {1, 2, 3, 4, 5};
    auto res = row_slope2(a.data(), a.data(), a.size());
    REQUIRE(res == Approx(1));
}

TEST_CASE("row_diff_mad", "[row_diff_mad]") {
    std::vector<double> a = {1, 2, 3, 4, 5};
    auto res = row_diff_mad(a.data(), a.size());
    REQUIRE(res == Approx(1));
}

TEST_CASE("row_sum_ad", "[row_sum_ad]") {
    std::vector<double> a = {1, 2, 3, 4, 5};
    auto res = row_sum_ad(a.data(), a.size(), 3);
    REQUIRE(res == Approx(6));
}

TEST_CASE("row_sum_ad2", "[row_sum_ad2]") {
    std::vector<double> a = {1, 2, 3, 4, 5};
    auto res = row_sum_ad2(a.data(), a.data(), a.size(), 3);
    REQUIRE(res == Approx(1.2));
}

TEST_CASE("row_sum_ed", "[row_sum_ed]") {
    std::vector<double> a = {1, 2, 3, 4, 5};
    auto res = row_sum_ed(a.data(), a.size(), 3);
    REQUIRE(res == Approx(2.0064294488));
}

TEST_CASE("row_sum_ed2", "[row_sum_ed2]") {
    std::vector<double> a = {1, 2, 3, 4, 5};
    auto res = row_sum_ed2(a.data(), a.data(), a.size(), 3);
    REQUIRE(res == Approx(0.4012858898));
}
