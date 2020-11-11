#include <math_dummy.h>
#include <algorithm>
#include <cmath>
#include <cstddef>

namespace ornate {
int dummy_add(int x, int y) { return x + y; }
int universal_answer() { return 42; }
double dummy_ts_cross(const double *x_, const double *y, std::size_t i, std::size_t n) {
    const std::size_t size = i + 1;
    if (size < n + 2) {
        return NAN;
    }
    double x0 = x_[i - n - 1];
    double y0 = y[i - n - 1];
    if (std::isnan(x0) || std::isnan(y0)) {
        return NAN;
    } else {
        double x1_ = x_[i - n];
        double y1 = y[i - n];
        if (std::isnan(x1_) || std::isnan(y1)) {
            return NAN;
        } else if (x0 < y0 && x1_ > y1) {
            return 1;
        } else if (x0 > y0 && x1_ < y1) {
            return -1;
        } else {
            return 0;
        }
    }
}

double dummy_ts_backward_cpn(const double *x, std::size_t i, std::size_t n, int sign) {
    const std::size_t size = i + 1;
    if (size < n) {
        return NAN;
    }
    const std::size_t N = (n == 0) ? size : std::min(size, n);
    double num = 0;
    for (std::size_t j = 0; j < N; j++) {
        if (std::isnan(x[i - j]) || (sign == 0 ? x[i] : sign) * x[i - j] <= 0) {
            break;
        } else {
            num++;
        }
    }
    return num;
}

double dummy_ts_acp(const double *x, int i, int n, int lag) {
    int bg = i - n + 1;
    int end = i + 1;

    bg = std::max(bg, 0);
    std::size_t valid_cnt = 0;
    std::size_t positive_cnt = 0;
    for (int j = bg; j < end - lag; j++) {
        const double vx = x[j];
        const double vx2 = x[j + lag];
        if (std::isnan(vx) or std::isnan(vx2)) continue;
        valid_cnt++;
        double res = vx * vx2;
        if (res > 0) {
            positive_cnt++;
        }
    }

    if (valid_cnt < 1) return NAN;
    return positive_cnt * 1.0 / valid_cnt;
}
}  // namespace ornate
