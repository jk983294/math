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
    std::normal_distribution<double> nd1(0., 0.008);
    std::normal_distribution<double> nd2(0., 0.004);

    int n = 100;
    int ins_num = 1;
    vector<double> ret(n * ins_num, 0);
    vector<double> signal(n * ins_num, 0);
    for (int i = 0; i < n * ins_num; ++i) {
        ret[i] = nd(generator);
        if (ret[i] > 0) {
            signal[i] = ret[i] + nd2(generator);
        } else {
            signal[i] = ret[i] + nd1(generator);
        }
    }

    double signal_median = ornate::median(signal);
    cout << "median=" << signal_median << endl;
    printf("original %f,%f,%f\n", ornate::rcor(ret, signal), ornate::rcor(ret, signal, 1), ornate::corr(ret, signal));

    auto signal_minus_median = signal;
    ornate::vs_minus_inplace(signal_minus_median, signal_median);
    printf("signal_minus_median %f,%f,%f\n", ornate::rcor(ret, signal_minus_median),
           ornate::rcor(ret, signal_minus_median, 1), ornate::corr(ret, signal_minus_median));

    auto opposite_signal = signal;
    ornate::vs_multiply_inplace(opposite_signal, -1.);
    ornate::vs_multiply_inplace(opposite_signal, ornate::median(opposite_signal));
    printf("opposite_signal %f,%f,%f\n", ornate::rcor(ret, opposite_signal), ornate::rcor(ret, opposite_signal, 1),
           ornate::corr(ret, opposite_signal));
    return 0;
}
