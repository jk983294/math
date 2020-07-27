#include <math_random.h>
#include <math_stats_rolling.h>
#include <cstdlib>
#include <iostream>

using namespace ornate;
using namespace std;

/**
 * for double 1e-14 count=195351, v1=0.357644, v2=0.357644
 * for float 1e-5 count=89575, v1=0.747054, v2=0.747044
 */

double mean(const std::deque<double>& m_data) {
    double sum = 0;
    for (auto v : m_data) {
        sum += v;
    }
    return sum / m_data.size();
}

using TestPrecision = double;

int main() {
    mean_rolling<TestPrecision> sr(4);

    size_t count = 0;
    while (true) {
        auto vec = generate_uniform_float<TestPrecision>(100000, 0.0, 1.0);
        for (TestPrecision v : vec) {
            TestPrecision accum_value = sr(v);
            TestPrecision exact_ans = mean(sr.m_data);
            if (fabs(accum_value - exact_ans) > 1e-14) {
                printf("count=%zu, v1=%f, v2=%f", count, accum_value, exact_ans);
                return 0;
            }
            ++count;
        }
    }
    return 0;
}
