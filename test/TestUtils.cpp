#include <math_utils.h>
#include <iostream>
#include "catch.hpp"
#include "math_stats.h"

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
