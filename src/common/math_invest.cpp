#include <math_indicator.h>
#include <vector>

namespace ornate {
std::vector<double> t0_nav(const double* signals, const double* ret1, size_t len,
                           double long_t, double short_t, double cost) {
    std::vector<double> navs(len, 1);
    if (len <= 1) return navs;

    double nav = 1.;
    for (size_t i = 0; i < len - 1; ++i) {
        if (std::isfinite(signals[i]) && std::isfinite(ret1[i])) {
            if (std::isfinite(long_t) && signals[i] >= long_t) {
                nav *= (1.0 + ret1[i] - cost);
            } else if (std::isfinite(short_t) && signals[i] <= short_t) {
                nav *= (1.0 - ret1[i] - cost);
            }
        }
        navs[i + 1] = nav;
    }
    return navs;
}
}