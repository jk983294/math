#include <math_random.h>
#include <math_stats.h>
#include <iostream>

using namespace ornate;
using namespace std;

int main() {
    //    std::seed_seq seed{42};
    //    default_random_engine generator{seed};
    random_device rd;
    mt19937 generator(rd());
    std::normal_distribution<double> nd(0., 0.01);

    int n = 100;
    int ins_num = 1;
    vector<double> ret(n * ins_num, 0);
    for (int i = 0; i < n * ins_num; ++i) {
        ret[i] = nd(generator);
    }

    ret[0] *= 100;
    ret[1] *= 10;
    ret[2] *= 5;

    double sd_, mean_, change_ratio;
    std::tie(sd_, mean_, change_ratio) = ornate::lr_reduce_outlier(ret);
    printf("%f,%f,%f\n", sd_, mean_, change_ratio);
    return 0;
}
