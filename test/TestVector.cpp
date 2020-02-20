#include <math_utils.h>
#include <iostream>
#include "catch.hpp"
#include "math_vector.h"

using namespace std;
using namespace ornate;

TEST_CASE("vs inplace", "[MathVector]") {
    vector<int> x{1, 2, 3};
    vs_multiply_inplace(x, 2);
    REQUIRE(x[0] == 2);
    REQUIRE(x[2] == 6);
}

TEST_CASE("vs non inplace", "[MathVector]") {
    vector<int> x{1, 2, 3};
    auto y = vs_multiply(x, 2);
    REQUIRE(y[0] == 2);
    REQUIRE(y[2] == 6);
}

TEST_CASE("vv non inplace", "[MathVector]") {
    vector<int> x{1, 2, 3};
    vector<int> y{1, 2, 3};
    auto z = vv_multiply(x, y);
    REQUIRE(z[0] == 1);
    REQUIRE(z[2] == 9);
}

TEST_CASE("non inplace", "[MathVector]") {
    vector<int> x = linspace(1, 5, 5);
    REQUIRE(x[0] == 1);
    REQUIRE(x[4] == 5);
}

TEST_CASE("get sorted index", "[MathVector]") {
    vector<int> x{3, 1, 2};
    auto ret = get_sorted_index(x);
    REQUIRE(ret[0] == 1);
    REQUIRE(ret[1] == 2);
    REQUIRE(ret[2] == 0);

    ret = get_sorted_index(x, false);
    REQUIRE(ret[0] == 0);
    REQUIRE(ret[1] == 2);
    REQUIRE(ret[2] == 1);
}

TEST_CASE("get sorted rank", "[MathVector]") {
    vector<int> x{3, 1, 2};
    auto ret = get_sorted_rank(x);
    REQUIRE(ret[0] == 2);
    REQUIRE(ret[1] == 0);
    REQUIRE(ret[2] == 1);

    ret = get_sorted_rank(x, false);
    REQUIRE(ret[0] == 0);
    REQUIRE(ret[1] == 2);
    REQUIRE(ret[2] == 1);
}

TEST_CASE("slice", "[MathVector]") {
    vector<int> x{1, 2, 3};
    auto ret = slice(x, 0, 1);
    REQUIRE(ret.size() == 1);
    REQUIRE(ret[0] == 1);

    ret = slice(x, 0, 2);
    REQUIRE(ret.size() == 2);
    REQUIRE(ret[0] == 1);
    REQUIRE(ret[1] == 2);

    ret = slice(x, 0, 0);
    REQUIRE(ret.size() == 3);

    ret = slice(x, -1, 0);
    REQUIRE(ret.size() == 1);
    REQUIRE(ret[0] == 3);
}

TEST_CASE("filter", "[MathVector]") {
    vector<float> x{1, 2, NAN, 3, NAN};
    vector<uint32_t> index(x.size());
    auto num = filter(x.data(), x.size(), index.data());
    REQUIRE(num == 3);
    REQUIRE(index[0] == 0);
    REQUIRE(index[1] == 1);
    REQUIRE(index[2] == 3);

    x = {1, 2, NAN, 3, NAN};
    vector<uint32_t> idx = filter(x);
    REQUIRE(idx.size() == 3);
    REQUIRE(idx[0] == 0);
    REQUIRE(idx[1] == 1);
    REQUIRE(idx[2] == 3);
}

TEST_CASE("rank", "[MathVector]") {
    vector<float> x{3, 2, NAN, 1, NAN};
    ornate::rank(x);
    REQUIRE(x.size() == 5);
    REQUIRE(x[0] == 1.0);
    REQUIRE(x[1] == 0.5);
    REQUIRE(x[3] == 0.0);
}

TEST_CASE("normalize", "[MathVector]") {
    vector<float> x{3, 2, NAN, 1, NAN};
    ornate::normalize(x);
    REQUIRE(x.size() == 5);
    REQUIRE(x[0] == 1.0);
    REQUIRE(x[1] == 0.0);
    REQUIRE(x[3] == -1.0);
}

TEST_CASE("powerf", "[MathVector]") {
    vector<float> x{3, 2, NAN, 1, NAN};
    ornate::powerf(x, 2.0f);
    REQUIRE(x.size() == 5);
    REQUIRE(x[0] == 9.0);
    REQUIRE(x[1] == 4.0);
    REQUIRE(x[3] == 1.0);
}

TEST_CASE("power rank_first", "[MathVector]") {
    vector<float> x{3, 2, NAN, 1, NAN};
    ornate::power(x, 2.0f, true);
    REQUIRE(x.size() == 5);
    REQUIRE(x[0] == 0.25);
    REQUIRE(x[1] == 0.0);
    REQUIRE(x[3] == -0.25);
}

TEST_CASE("sum", "[MathVector]") {
    vector<float> x{3, 2, NAN, 1, NAN};
    REQUIRE(ornate::sum(x) == 6.0);
}
