#pragma once

#include <cmath>
#include <cstddef>

namespace ornate {
double row_sum(const double* a, int len);
double row_sum2(const double* a, const double* b, int len);
int row_cut(const double* a, double threshold, int len);
int row_sum_cut(const double* a, double threshold, int len);
double row_mean(const double* a, int len);
double row_mean2(const double* a, const double* b, int len);
double row_sd(const double* a, int len);
double row_sd2(const double* a, const double* b, int len);
double row_ema(const double* a, int len, double hl);
double row_ema2(const double* a, const double* b, int len, double hl);
double row_median(const double* a, int len);
double row_slope(const double* a, int len, bool xnorm = false, bool ynorm = false, bool intercept = false);
double row_slope2(const double* a, const double* b, int len, bool xnorm = false, bool ynorm = false,
                  bool intercept = false);
double row_diff_mad(const double* a, int len);
double row_sum_ad(const double* a, int len, double center);
double row_sum_ad2(const double* a, const double* b, int len, double center);
double row_sum_ed(const double* a, int len, double center, double lambda = 1.0);
double row_sum_ed2(const double* a, const double* b, int len, double center, double lambda = 1.0);
double row_max(const double* a, int len);
double row_min(const double* a, int len);
int row_which_max(const double* a, int len);
int row_which_min(const double* a, int len);
int row_which_eq(const double* a, double val, int len);
int row_which_gt(const double* a, double val, int len);
int row_which_lt(const double* a, double val, int len);
int row_which_gte(const double* a, double val, int len);
int row_which_lte(const double* a, double val, int len);
double row_extract(const double* a, int idx, int len);
}  // namespace ornate