#include <math_segment.h>
#include "catch.hpp"

using namespace std;
using namespace ornate;

TEST_CASE("segment_extract_skip_ti", "[segment]") {
    std::vector<int> data = {0, 1, 2, 3, 4, 5};
    auto skiped = ornate::segment_extract_skip_ti(data.data(), data.size(), 2, 3, {1});
    REQUIRE(skiped.size() == 4);
    REQUIRE(skiped.front() == 0);
    REQUIRE(skiped.back() == 5);

    skiped = ornate::segment_extract_skip_ti(data.data(), data.size(), 2, 3, {0, 2});
    REQUIRE(skiped.size() == 2);
    REQUIRE(skiped.front() == 2);
    REQUIRE(skiped.back() == 3);
}
