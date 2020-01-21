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
