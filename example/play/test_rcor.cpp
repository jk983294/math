#include <math_random.h>
#include <math_stats.h>
#include <iostream>

using namespace ornate;
using namespace std;

static void demean_test();

int main() {
    //    std::seed_seq seed{42};
    //    default_random_engine generator{seed};
    random_device rd;
    mt19937 generator(rd());
    std::normal_distribution<double> nd(0., 0.01);
    std::normal_distribution<double> nd1(0., 0.004);

    int n = 100;
    int ins_num = 1;
    vector<double> ret(n * ins_num, 0);
    vector<double> signal(n * ins_num, 0);
    for (int i = 0; i < n * ins_num; ++i) {
        ret[i] = nd(generator);
        if (ret[i] > 0) {
            signal[i] = ret[i] + nd1(generator) + 1.1;
        } else {
            signal[i] = ret[i] + nd1(generator);
        }
    }

    auto opposite_signal = signal;
    ornate::vs_multiply_inplace(opposite_signal, -1.);
    double signal_median = ornate::median(signal);
    double signal_median2 = ornate::median(opposite_signal);
    cout << "median=" << signal_median << "," << signal_median2 << endl;
    cout << "skew=" << ornate::math_skew(signal) << endl;
    printf("%f,%f,%f\n", signal[0], signal[1], signal[2]);
    printf("original %f,%f,%f,%f\n", ornate::rcor(ret, signal), ornate::rcor(ret, signal, 1),
           ornate::rcor(ret, signal, -1), ornate::corr(ret, signal));

    auto s_sub_median = signal;
    ornate::vs_minus_inplace(s_sub_median, signal_median);
    printf("%f,%f,%f\n", s_sub_median[0], s_sub_median[1], s_sub_median[2]);
    printf("signal_minus_median %f,%f,%f,%f\n", ornate::rcor(ret, s_sub_median), ornate::rcor(ret, s_sub_median, 1),
           ornate::rcor(ret, s_sub_median, -1), ornate::corr(ret, s_sub_median));

    ornate::vs_minus_inplace(opposite_signal, signal_median2);
    printf("%f,%f,%f\n", opposite_signal[0], opposite_signal[1], opposite_signal[2]);
    printf("opposite_signal %f,%f,%f,%f\n", ornate::rcor(ret, opposite_signal), ornate::rcor(ret, opposite_signal, 1),
           ornate::rcor(ret, opposite_signal, -1), ornate::corr(ret, opposite_signal));

    demean_test();
    return 0;
}

/**
 * pcor can handle linear drift, it demean automatically
 * rcor needs signal has good shape, can not drift, it doesn't demean
 */
void demean_test() {
    random_device rd;
    mt19937 generator(rd());
    std::normal_distribution<double> nd(0., 0.01);
    std::normal_distribution<double> nd1(0., 0.004);

    int n = 100;
    int ins_num = 1;
    vector<double> ret(n * ins_num, 0);
    vector<double> signal(n * ins_num, 0);
    vector<double> signal1(n * ins_num, 0);
    for (int i = 0; i < n * ins_num; ++i) {
        ret[i] = nd(generator);
        signal1[i] = ret[i] + nd1(generator);
        signal[i] = signal1[i] + 1.1;
    }

    printf("signal_minus_median %f,%f,%f,%f\n", ornate::rcor(ret, signal), ornate::rcor(ret, signal1),
           ornate::corr(ret, signal), ornate::corr(ret, signal1));
}
