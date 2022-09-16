#ifndef ORNATE_MATH_DUMMY_H
#define ORNATE_MATH_DUMMY_H

#include <cmath>
#include <cstddef>
#include <functional>
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

namespace dminner {
double mean(const double *x, int start, int end, double fill);

double quantile(const double *x, double q, int start, int end, double fill);

double max(const double *x, int start, int end, double fill);

double min(const double *x, int start, int end, double fill);
}

using COMPARE = std::function<bool(double, double)>;
using COND_FUNC = std::function<double(COMPARE, int, const double *, const double *, double, int, int, double)>;

double cond_mean(const COMPARE& g, int method, const double *x, const double *y,
                 double q, int start, int end, double fill);

double cond_max(const COMPARE& g, int method, const double *x, const double *y,
                double q, int start, int end, double fill);

double cond_min(const COMPARE& g, int method, const double *x, const double *y,
                double q, int start, int end, double fill);

double cond_sd(const COMPARE& g, int method, const double *x, const double *y,
               double q, int start, int end, double fill);

std::vector<double> ts_gte_mean(const std::vector<double> &cond, const std::vector<double> &value,
                                int n, double q = 0.5, double fill = NAN, int method = 1,
                                int least = 3, bool partial = false);

std::vector<double> ts_lte_mean(const std::vector<double> &cond, const std::vector<double> &value,
                                int n, double q = 0.5, double fill = NAN, int method = 1,
                                int least = 3, bool partial = false);

std::vector<double> ts_gte_max(const std::vector<double> &cond, const std::vector<double> &value,
                               int n, double q = 0.5, double fill = NAN, int method = 1,
                               int least = 3, bool partial = false);

std::vector<double> ts_lte_max(const std::vector<double> &cond, const std::vector<double> &value,
                               int n, double q = 0.5, double fill = NAN, int method = 1,
                               int least = 3, bool partial = false);

std::vector<double> ts_gte_min(const std::vector<double> &cond, const std::vector<double> &value,
                               int n, double q = 0.5, double fill = NAN, int method = 1,
                               int least = 3, bool partial = false);

std::vector<double> ts_lte_min(const std::vector<double> &cond, const std::vector<double> &value,
                               int n, double q = 0.5, double fill = NAN, int method = 1,
                               int least = 3, bool partial = false);

std::vector<double> ts_gte_sd(const std::vector<double> &cond, const std::vector<double> &value,
                              int n, double q = 0.5, double fill = NAN, int method = 1,
                              int least = 3, bool partial = false);

std::vector<double> ts_lte_sd(const std::vector<double> &cond, const std::vector<double> &value,
                              int n, double q = 0.5, double fill = NAN, int method = 1,
                              int least = 3, bool partial = false);

std::vector<double> ts_gte_sum(const std::vector<double> &cond, const std::vector<double> &value,
                               int n, double q = 0.5, double fill = NAN, int method = 1,
                               int least = 3, bool partial = false);

std::vector<double> ts_lte_sum(const std::vector<double> &cond, const std::vector<double> &value,
                               int n, double q = 0.5, double fill = NAN, int method = 1,
                               int least = 3, bool partial = false);
}  // namespace ornate

#endif
