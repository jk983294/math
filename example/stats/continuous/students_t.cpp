#include <zerg_stl.h>
#include <iostream>
#include <stats.hpp>
#include <vector>

using namespace std;

/**
 * Student's t distribution, continuous
 * https://en.wikipedia.org/wiki/Student%27s_t-distribution
 *
 * f(x;ν)
 */

int main() {
    double sample_count = 14.0;

    /**
     * f(x;ν)
     */
    cout << "density " << stats::dt(1.0, sample_count) << endl;  // 0.233581
    cout << "density " << stats::dt(0.5, sample_count) << endl;  // 0.343171

    /**
     * P(X <= x)
     */
    cout << "cdf " << stats::pt(1.0, sample_count) << endl;  // 0.832859
    cout << "cdf " << stats::pt(0, sample_count) << endl;    // 0.5

    /**
     * find x such that cdf(x) = a
     */
    cout << "quantile " << stats::qt(0.841345, sample_count) << endl;  // 1.03701
    cout << "quantile " << stats::qt(0.5, sample_count) << endl;       // 0
    cout << "quantile " << stats::qt(0.4, sample_count) << endl;       // -0.258213
    cout << "quantile " << stats::qt(0.3, sample_count) << endl;       // -0.536552
}
