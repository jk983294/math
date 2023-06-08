#include <math_random.h>
#include <math_stats_once.h>
#include <chrono>
#include <iostream>

using namespace ornate;
using namespace std;
using namespace std::chrono;

int main() {
    //    std::seed_seq seed{42};
    //    default_random_engine generator{seed};
    random_device rd;
    mt19937 generator(rd());
    std::normal_distribution<double> nd(0., 1);

    QuantileOnce qo(0.01);

    int n = 10000000;
    vector<double> ret(n, 0);
    for (int i = 0; i < n; ++i) {
        ret[i] = nd(generator);
    }

    steady_clock::time_point t1 = steady_clock::now();
    for (int i = 0; i < n; ++i) {
        qo.insert(ret[i]);
    }
    printf("%f,%f,%f,%f,%f,%f,%f,%f,%f\n", qo.query(0.1), qo.query(0.2), qo.query(0.3), qo.query(0.4), qo.query(0.5),
           qo.query(0.6), qo.query(0.7), qo.query(0.8), qo.query(0.9));
    steady_clock::time_point t2 = steady_clock::now();

    std::sort(ret.begin(), ret.end());
    printf("%f,%f,%f,%f,%f,%f,%f,%f,%f\n", ret[0.1 * n], ret[0.2 * n], ret[0.3 * n], ret[0.4 * n], ret[0.5 * n],
           ret[0.6 * n], ret[0.7 * n], ret[0.8 * n], ret[0.9 * n]);
    steady_clock::time_point t3 = steady_clock::now();

    printf("rank:\n");
    for (int j = 0; j < 10; ++j) {
        double val = ret[(0.1 * j + 0.05) * n];
        printf("%f,%f\n", val, qo.query_rank(val));
    }

    long tc = nanoseconds{t2 - t1}.count();
    long tc1 = nanoseconds{t3 - t2}.count();
    printf("%ld,%ld\n", tc, tc1);

    return 0;
}
