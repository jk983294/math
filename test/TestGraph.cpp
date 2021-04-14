#include <math_utils.h>
#include <iostream>
#include "catch.hpp"
#include "math_graph.h"

using namespace std;
using namespace ornate;

TEST_CASE("splitter 1", "[GraphSplitter]") {
    GraphSplitter splitter(9, 4);
    splitter.add_constrain({0, 1, 5});
    splitter.add_constrain({0, 1, 6});
    splitter.add_constrain({0, 1, 3});
    splitter.add_constrain({0, 1, 4});
    splitter.add_constrain({0, 2, 7});
    splitter.add_constrain({0, 2, 8});
    splitter.add_constrain({0, 2, 3});
    splitter.add_constrain({0, 2, 4});
    splitter.add_constrain({1, 2, 5});
    splitter.add_constrain({1, 2, 6});
    splitter.add_constrain({1, 2, 7});
    splitter.add_constrain({1, 2, 8});

    bool status = splitter.work();
    REQUIRE(status);
    REQUIRE(splitter.m_result_constrains.size() == 6);

    splitter.m_max_group_nodes = 5;
    status = splitter.work();
    REQUIRE(status);
    REQUIRE(splitter.m_result_constrains.size() == 4);

    splitter.m_max_group_nodes = 6;
    status = splitter.work();
    REQUIRE(status);
    REQUIRE(splitter.m_result_constrains.size() == 2);

    splitter.m_max_group_nodes = 9;
    status = splitter.work();
    REQUIRE(status);
    REQUIRE(splitter.m_result_constrains.size() == 1);

    splitter.m_max_group_nodes = 3;
    status = splitter.work();
    REQUIRE(status);
    REQUIRE(splitter.m_result_constrains.size() == 12);

    splitter.m_max_group_nodes = 2;
    status = splitter.work();
    REQUIRE(!status);
    REQUIRE(splitter.m_result_constrains.size() == 0);
}

TEST_CASE("splitter 2", "[GraphSplitter]") {
    GraphSplitter splitter(9, 1);
    splitter.add_constrain({0, 1});
    splitter.add_constrain({0, 2});
    splitter.add_constrain({1, 2});

    splitter.add_constrain({0, 3});
    splitter.add_constrain({0, 4});
    splitter.add_constrain({1, 5});
    splitter.add_constrain({1, 6});
    splitter.add_constrain({2, 7});
    splitter.add_constrain({2, 8});

    bool status = splitter.work();
    REQUIRE(!status);
    REQUIRE(splitter.m_result_constrains.size() == 0);

    splitter.m_max_group_nodes = 2;
    status = splitter.work();
    REQUIRE(status);
    REQUIRE(splitter.m_result_constrains.size() == 9);

    splitter.m_max_group_nodes = 3;
    status = splitter.work();
    REQUIRE(status);
    REQUIRE(splitter.m_result_constrains.size() == 4);

    splitter.m_max_group_nodes = 4;
    status = splitter.work();
    REQUIRE(status);
    REQUIRE(splitter.m_result_constrains.size() == 3);

    splitter.m_max_group_nodes = 5;
    status = splitter.work();
    REQUIRE(status);
    REQUIRE(splitter.m_result_constrains.size() == 3);

    splitter.m_max_group_nodes = 6;
    status = splitter.work();
    REQUIRE(status);
    REQUIRE(splitter.m_result_constrains.size() == 2);

    splitter.m_max_group_nodes = 8;
    status = splitter.work();
    REQUIRE(status);
    REQUIRE(splitter.m_result_constrains.size() == 2);

    splitter.m_max_group_nodes = 9;
    status = splitter.work();
    REQUIRE(status);
    REQUIRE(splitter.m_result_constrains.size() == 1);

    splitter.m_max_group_nodes = 10;
    status = splitter.work();
    REQUIRE(status);
    REQUIRE(splitter.m_result_constrains.size() == 1);
}

TEST_CASE("splitter 3", "[GraphSplitter]") {
    GraphSplitter splitter(9, 1);
    splitter.add_constrain({0, 3});
    splitter.add_constrain({0, 4});
    splitter.add_constrain({1, 5});
    splitter.add_constrain({1, 6});
    splitter.add_constrain({2, 7});
    splitter.add_constrain({2, 8});

    bool status = splitter.work();
    REQUIRE(!status);
    REQUIRE(splitter.m_result_constrains.size() == 0);

    splitter.m_max_group_nodes = 2;
    status = splitter.work();
    REQUIRE(status);
    REQUIRE(splitter.m_result_constrains.size() == 6);

    splitter.m_max_group_nodes = 3;
    status = splitter.work();
    REQUIRE(status);
    REQUIRE(splitter.m_result_constrains.size() == 3);

    splitter.m_max_group_nodes = 4;
    status = splitter.work();
    REQUIRE(status);
    REQUIRE(splitter.m_result_constrains.size() == 3);

    splitter.m_max_group_nodes = 5;
    status = splitter.work();
    REQUIRE(status);
    REQUIRE(splitter.m_result_constrains.size() == 2);

    splitter.m_max_group_nodes = 6;
    status = splitter.work();
    REQUIRE(status);
    REQUIRE(splitter.m_result_constrains.size() == 2);

    splitter.m_max_group_nodes = 8;
    status = splitter.work();
    REQUIRE(status);
    REQUIRE(splitter.m_result_constrains.size() == 2);

    splitter.m_max_group_nodes = 9;
    status = splitter.work();
    REQUIRE(status);
    REQUIRE(splitter.m_result_constrains.size() == 1);

    splitter.m_max_group_nodes = 10;
    status = splitter.work();
    REQUIRE(status);
    REQUIRE(splitter.m_result_constrains.size() == 1);
}
