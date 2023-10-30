#include <math_indicator.h>
#include <vector>
#include <algorithm>

namespace ornate {
double HHS(const double* close, const double* volume, const double* oi, std::size_t i, std::size_t n) {
    std::size_t s_pos = 0;
    if (i >= n) s_pos = i - n;
    std::vector<double> vals;
    vals.reserve(n);
    for (std::size_t j = s_pos; j <= i; ++j) {
        if (std::isfinite(volume[j])) {
            vals.push_back(volume[j]);
        }
    }

    if (vals.empty()) return NAN;
    std::sort(vals.begin(), vals.end());
    int valid_n = vals.size();
    double one_third =  vals.front() + (vals.back() - vals.front()) / 3.;
    double two_third =  vals.back() - (vals.back() - vals.front()) / 3.;
    if (valid_n > 1 && valid_n <= 3) {
        one_third = vals.front() + (vals[i] - vals.front()) / 2.;
        two_third = vals.back() - (vals.back() - *(vals.end() - 2)) / 2.;
    } else if (valid_n > 3) {
        one_third = vals[std::lrint(std::ceil(valid_n / 3.))];
        two_third = vals[std::lrint(std::floor(2 * valid_n / 3.))];
    }

    int bull = 0, bear = 0;
    for (std::size_t j = s_pos; j <= i; ++j) {
        if (not std::isfinite(close[j]) || j == 0 || not std::isfinite(close[j - 1])) continue;
        if (not std::isfinite(volume[j]) || not std::isfinite(volume[j - 1])) continue;
        if (close[j] > close[j - 1]) {
            if (volume[j] >= two_third) bull += 1;
            else if (volume[j] <= one_third) bear += 1;

            if (oi) {
                if (oi[j] > oi[j - 1]) bull += 1;
                else bear += 1;
            }
        } else if (close[j] < close[j - 1]) {
            if (volume[j] >= two_third) bear += 1;
            else if (volume[j] <= one_third) bull += 1;

            if (oi) {
                if (oi[j] > oi[j - 1]) bear += 1;
                else bull += 1;
            }
        }
    }
    return bull - bear;
}
}