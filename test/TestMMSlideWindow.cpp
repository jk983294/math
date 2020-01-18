#include <math_utils.h>
#include <iostream>
#include "MMSlideWindow.h"
#include "catch.hpp"

using namespace std;
using namespace ornate;

TEST_CASE("window", "[MMSlideWindow]") {
    REQUIRE(ornate::IsValidData(1));
    REQUIRE(!ornate::IsValidData(NAN));
    REQUIRE(!ornate::IsValidData(INFINITY));

    MaxMinSlideWindow<int> window(3);
    window.set_nan(-42);
    window.add_new(3);
    REQUIRE(window.min_val() == 3);
    REQUIRE(window.max_val() == 3);

    window.add_new(5);
    window.print_all_data();
    REQUIRE(window.min_val() == 3);
    REQUIRE(window.max_val() == 5);

    window.add_new(1);
    window.print_all_data();
    REQUIRE(window.min_val() == 1);
    REQUIRE(window.max_val() == 5);

    REQUIRE(window.mink_val(0) == 1);
    REQUIRE(window.mink_val(1) == 3);
    REQUIRE(window.mink_val(2) == 5);

    REQUIRE(window.maxk_val(0) == 5);
    REQUIRE(window.maxk_val(1) == 3);
    REQUIRE(window.maxk_val(2) == 1);

    REQUIRE(window.mink_pos(0) == 2);
    REQUIRE(window.mink_pos(1) == 0);
    REQUIRE(window.mink_pos(2) == 1);

    REQUIRE(window.maxk_pos(0) == 1);
    REQUIRE(window.maxk_pos(1) == 0);
    REQUIRE(window.maxk_pos(2) == 2);

    window.add_new(2);
    REQUIRE(window.min_val() == 1);
    REQUIRE(window.max_val() == 5);

    window.add_new(4);
    REQUIRE(window.min_val() == 1);
    REQUIRE(window.max_val() == 4);

    REQUIRE(window.mink_val(0) == 1);
    REQUIRE(window.mink_val(1) == 2);
    REQUIRE(window.mink_val(2) == 4);

    REQUIRE(window.maxk_val(0) == 4);
    REQUIRE(window.maxk_val(1) == 2);
    REQUIRE(window.maxk_val(2) == 1);

    REQUIRE(window.mink_pos(0) == 0);
    REQUIRE(window.mink_pos(1) == 1);
    REQUIRE(window.mink_pos(2) == 2);

    REQUIRE(window.maxk_pos(0) == 2);
    REQUIRE(window.maxk_pos(1) == 1);
    REQUIRE(window.maxk_pos(2) == 0);

    REQUIRE(window.max_pos() == 2);
    REQUIRE(window.min_pos() == 0);
    REQUIRE(window.size() == 3);

    window.clear();
    REQUIRE(window.size() == 0);
}
