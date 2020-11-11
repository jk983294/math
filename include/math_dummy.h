#ifndef ORNATE_MATH_DUMMY_H
#define ORNATE_MATH_DUMMY_H

#include <cstddef>

namespace ornate {
int dummy_add(int x, int y);
int universal_answer();
int dummy_prod(int x, int y);
double dummy_ts_cross(const double *x_, const double *y, std::size_t i, std::size_t n);
double dummy_ts_backward_cpn(const double *x, std::size_t i, std::size_t n, int sign);
double dummy_ts_acp(const double *x, int i, int n, int lag);
}  // namespace ornate

#endif
