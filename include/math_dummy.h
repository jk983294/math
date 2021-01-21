#ifndef ORNATE_MATH_DUMMY_H
#define ORNATE_MATH_DUMMY_H

#include <cmath>
#include <cstddef>
#include <vector>

namespace ornate {
int dummy_add(int x, int y);
int universal_answer();
int dummy_prod(int x, int y);
double dummy_ts_cross(const double* x_, const double* y, std::size_t i, std::size_t n);
double dummy_ts_argmax(const double* x_, std::size_t i, std::size_t n);
double dummy_ts_backward_cpn(const double* x, std::size_t i, std::size_t n, int sign);
double dummy_ts_acp(const double* x, int i, int n, int lag);
double dummy_dcor(const std::vector<double>& x_, const std::vector<double>& y_);
double dummy_r2(const std::vector<double>& x_, const std::vector<double>& y_, double a, double b, int window);
double dummy_r2(const std::vector<double>& x1_, const std::vector<double>& x2_, const std::vector<double>& y_,
                double b0, double b1, double b2, int window);
double dummy_r2_no_slope(const std::vector<double>& x_, const std::vector<double>& y_, double b, int window);
double dummy_r2_no_slope(const std::vector<double>& x1_, const std::vector<double>& x2_, const std::vector<double>& y_,
                         double b1, double b2, int window);
double dummy_skew(const double* data_, int num);
double dummy_kurtosis(const double* data_, int num);
double dummy_skew(const std::vector<double>& data_);
double dummy_kurtosis(const std::vector<double>& data_);

template <typename T>
double dummy_decay(const std::vector<T>& data_) {
    int size = (int)data_.size();
    double res = 0;
    int count = 0;
    for (int i = 0; i < size; ++i) {
        if (std::isfinite(data_[i])) {
            res += data_[i] * (i + 1);
            count += (i + 1);
        }
    }
    if (count > 0) {
        return res / count;
    } else {
        return NAN;
    }
}
}  // namespace ornate

#endif
