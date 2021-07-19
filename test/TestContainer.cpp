#include <math_container.h>
#include <numeric>
#include "catch.hpp"

using namespace std;
using namespace ornate;

TEST_CASE("container2 data", "[container]") {
    rolling_data_container2<double> c2(2, 3, 1);
    vector<double> datum(1, 0);
    c2.push(datum);
    REQUIRE(c2.get_old_row0() == nullptr);
    REQUIRE(c2.get_old_row1() == nullptr);
    REQUIRE(*c2.get_new_row() == 0);

    datum.front() = 1;
    c2.push(datum);
    REQUIRE(c2.get_old_row0() == nullptr);
    REQUIRE(c2.get_old_row1() == nullptr);
    REQUIRE(*c2.get_new_row() == 1);

    datum.front() = 2;
    c2.push(datum);
    REQUIRE(*c2.get_old_row0() == 0);
    REQUIRE(c2.get_old_row1() == nullptr);
    REQUIRE(*c2.get_new_row() == 2);

    datum.front() = 3;
    c2.push(datum);
    REQUIRE(*c2.get_old_row0() == 1);
    REQUIRE(*c2.get_old_row1() == 0);
    REQUIRE(*c2.get_new_row() == 3);
}

TEST_CASE("container2 pointer", "[container]") {
    rolling_pointer_container2<double> c2(2, 3, 1);
    vector<double> datum(5);
    std::iota(datum.begin(), datum.end(), 0);
    c2.push(datum.data());
    REQUIRE(c2.get_old_row0() == nullptr);
    REQUIRE(c2.get_old_row1() == nullptr);
    REQUIRE(*c2.get_new_row() == 0);

    c2.push(datum.data() + 1);
    REQUIRE(c2.get_old_row0() == nullptr);
    REQUIRE(c2.get_old_row1() == nullptr);
    REQUIRE(*c2.get_new_row() == 1);

    c2.push(datum.data() + 2);
    REQUIRE(*c2.get_old_row0() == 0);
    REQUIRE(c2.get_old_row1() == nullptr);
    REQUIRE(*c2.get_new_row() == 2);

    c2.push(datum.data() + 3);
    REQUIRE(*c2.get_old_row0() == 1);
    REQUIRE(*c2.get_old_row1() == 0);
    REQUIRE(*c2.get_new_row() == 3);
}
