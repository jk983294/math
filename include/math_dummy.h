#ifndef ORNATE_MATH_DUMMY_H
#define ORNATE_MATH_DUMMY_H

#include <cstddef>
#include <vector>

namespace ornate {
int dummy_add(int x, int y);
int universal_answer();
int dummy_prod(int x, int y);
double dummy_ts_cross(const double* x_, const double* y, std::size_t i, std::size_t n);
double dummy_ts_backward_cpn(const double* x, std::size_t i, std::size_t n, int sign);
double dummy_ts_acp(const double* x, int i, int n, int lag);
double dummy_dcor(const std::vector<double>& x_, const std::vector<double>& y_);
double dummy_r2(const std::vector<double>& x_, const std::vector<double>& y_, double a, double b, int window);
double dummy_r2(const std::vector<double>& x1_, const std::vector<double>& x2_, const std::vector<double>& y_,
                double b0, double b1, double b2, int window);
double dummy_r2_no_slope(const std::vector<double>& x_, const std::vector<double>& y_, double b, int window);
double dummy_r2_no_slope(const std::vector<double>& x1_, const std::vector<double>& x2_, const std::vector<double>& y_,
                         double b1, double b2, int window);
}  // namespace ornate

#endif
