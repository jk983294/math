#include <math_random.h>
#include <iostream>
#include "catch.hpp"

using namespace std;
using namespace ornate;

TEST_CASE("choice", "[math_random]") {
    vector<int> data = {1, 3, 2, 4};
    vector<double> weight_accum = build_choice_vector(data);
    REQUIRE(weight_accum.size() == 4);
    REQUIRE(weight_accum[0] == 0.1);
    REQUIRE(weight_accum[1] == 0.4);
    REQUIRE(weight_accum[2] == 0.6);
    REQUIRE(weight_accum[3] == 1.0);

    REQUIRE(choice_with_accum_weight(weight_accum, 0.05) == 0);
    REQUIRE(choice_with_accum_weight(weight_accum, 0.1) == 1);
    REQUIRE(choice_with_accum_weight(weight_accum, 0.3) == 1);
    REQUIRE(choice_with_accum_weight(weight_accum, 0.4) == 2);
    REQUIRE(choice_with_accum_weight(weight_accum, 0.5) == 2);
    REQUIRE(choice_with_accum_weight(weight_accum, 0.6) == 3);
    REQUIRE(choice_with_accum_weight(weight_accum, 0.7) == 3);
    REQUIRE(choice_with_accum_weight(weight_accum, 1.0) == 3);
    REQUIRE(choice_with_accum_weight(weight_accum, 1.1) == 3);
}
