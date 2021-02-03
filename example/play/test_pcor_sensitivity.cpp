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
    vector<double> signal(n * ins_num, 0);
    for (int i = 0; i < n * ins_num; ++i) {
        ret[i] = nd(generator);
        signal[i] = ret[i];
    }

    /**
     * 完全一样的序列，将其中一个值变大，相关性迅速下降
     */
    cout << "0," << ornate::corr(ret, signal) << endl;
    for (int j = 0; j < 10; ++j) {
        signal[0] *= 2;
        cout << (j + 1) << "," << ornate::corr(ret, signal) << endl;
    }

    /**
     * 完全不同的序列，初始相关性很低，将相同位置的值变相同，再让极端值同步变大，相关性迅速上升
     * 这说明只要预测准大变化，相关性就变大
     */
    for (int i = 0; i < n * ins_num; ++i) {
        signal[i] = nd(generator);
    }
    signal[0] = ret[0];
    cout << "0," << ornate::corr(ret, signal) << endl;
    for (int j = 0; j < 10; ++j) {
        ret[0] *= 2;
        signal[0] = ret[0];
        cout << (j + 1) << "," << ornate::corr(ret, signal) << endl;
    }

    /**
     * 离群值越来越大，则skew越来越大，说明大部分值都小于mean，kurtosis越来越大，重尾, 极端值比较多
     * skew < 0, 大部分 > mean
     * skew > 0, 大部分 < mean
     * kurtosis > 0, 重尾, 极端值比较多
     * kurtosis < 0, 轻尾, 极端值比较少, 大部分靠中间
     */
    ret[0] = abs(nd(generator));
    double mean_ = ornate::mean(ret);
    cout << "0," << mean_ << "," << ornate::math_skew(ret) << "," << ornate::math_kurtosis(ret) << endl;
    for (int j = 0; j < 10; ++j) {
        ret[0] *= 2;
        mean_ = ornate::mean(ret);
        cout << (j + 1) << "," << mean_ << "," << ornate::math_skew(ret) << "," << ornate::math_kurtosis(ret) << endl;
    }
    return 0;
}
