#include <math_utils.h>
#include "catch.hpp"

using namespace std;
using namespace ornate;

TEST_CASE("valid", "[MathStats]") {
    REQUIRE(ornate::IsValidData(1));
    REQUIRE(!ornate::IsValidData(NAN));
    REQUIRE(!ornate::IsValidData(INFINITY));
}

TEST_CASE("equal", "[MathStats]") {
    float a = 16777216.0f;
    float b = 16777256.0f;
    double c = 0.1000000000000000055511151231257827021181583404541015625;
    double d = 0.1000000000000000195511151231257827021181583404541015625;
    REQUIRE(!FloatEqual(a, b));
    REQUIRE(FloatEqual(c, d));
    REQUIRE(FloatEqual(NAN, NAN));
    REQUIRE(FloatEqual(INFINITY, INFINITY));
}

TEST_CASE("div_down", "[MathStats]") {
    REQUIRE(FloatEqual(div_down(1234, 0.023, 1), 0.001234));
    REQUIRE(FloatEqual(div_down(2, 0.023, 1), 0.002000));
    REQUIRE(FloatEqual(div_down(0.2, 0.023, 1), 0.000200));
    REQUIRE(FloatEqual(div_down(0.02, 0.023, 1), 0.000200));
    REQUIRE(FloatEqual(div_down(0.002, 0.023, 1), 0.00200));
}
