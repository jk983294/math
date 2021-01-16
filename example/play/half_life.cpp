#include <math_random.h>
#include <math_stats.h>
#include <zerg_template.h>
#include <cstdlib>
#include <iostream>

using namespace ornate;
using namespace std;

int main() {
    //    std::seed_seq seed{42};
    //    default_random_engine generator{seed};
    random_device rd;
    mt19937 generator(rd());
    std::normal_distribution<double> nd(0., 0.001);

    int n = 10000;
    int ins_num = 2;
    int window = 60;
    vector<double> ret(n * ins_num, 0);
    vector<double> fwd_sum_ret(n * ins_num, 0);
    for (int i = 0; i < n * ins_num; ++i) {
        ret[i] = nd(generator);
    }
    for (int i = 0; i < n * ins_num; ++i) {
        for (int j = 1; j <= window; ++j) {
            if (i + j * ins_num < n * ins_num) {
                fwd_sum_ret[i] += ret[i + j * ins_num];
            } else {
                break;
            }
        }
    }

    printf("%s\n", ztool::head(ret, 12).c_str());
    printf("%s\n", ztool::head(fwd_sum_ret, 12).c_str());

    for (int i = 0; i < window; ++i) {
        double pcor = ornate::corr(fwd_sum_ret.data(), fwd_sum_ret.data() + i * ins_num, (n - i) * ins_num);
        printf("shift=%d, %f %f, pcor=%f\n", i, fwd_sum_ret[0], fwd_sum_ret[i * ins_num], pcor);
    }

    return 0;
}
