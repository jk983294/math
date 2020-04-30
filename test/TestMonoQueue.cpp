#include <math_utils.h>
#include "math_monoqueue.h"
#include "catch.hpp"
#include <queue>

using namespace std;
using namespace ornate;

TEST_CASE("monoqueue min", "[monoqueue min]") {
    MonotoneQueue<int> q(5);
    REQUIRE(q.Top() == int32_nan);
    REQUIRE(q.TopIndex() == -1);

    for (int i = 0; i < 5; ++i) {
        q.Push(i);
        REQUIRE(q.Top() == 0);
        REQUIRE(q.TopIndex() == i);
    }

    for (int i = 5; i < 10; ++i) {
        q.Push(i);
        REQUIRE(q.Top() == i - 4);
        REQUIRE(q.TopIndex() == 4);
    }

    q.Push(0);
    REQUIRE(q.Top() == 0);
    REQUIRE(q.TopIndex() == 0);

    for (int i = 1; i < 5; ++i) {
        q.Push(i);
        REQUIRE(q.Top() == 0);
        REQUIRE(q.TopIndex() == i);
    }

    q.Push(-1);
    REQUIRE(q.Top() == -1);
    REQUIRE(q.TopIndex() == 0);
}

TEST_CASE("monoqueue max", "[monoqueue max]") {
    MonotoneQueue<int, std::less<>> q(5);
    REQUIRE(q.Top() == int32_nan);
    REQUIRE(q.TopIndex() == -1);

    for (int i = 4; i >= 0; --i) {
        q.Push(i);
        REQUIRE(q.Top() == 4);
        REQUIRE(q.TopIndex() == 4 - i);
    }

    for (int i = -1; i >= -5; --i) {
        q.Push(i);
        REQUIRE(q.Top() == i + 4);
        REQUIRE(q.TopIndex() == 4);
    }

    q.Push(0);
    REQUIRE(q.Top() == 0);
    REQUIRE(q.TopIndex() == 0);
}
