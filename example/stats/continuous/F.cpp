#include <zerg_stl.h>
#include <iostream>
#include <stats.hpp>
#include <vector>

using namespace std;

/**
 * F distribution, continuous
 * https://en.wikipedia.org/wiki/F-distribution
 * A random variate of the F-distribution with parameters d1 and d2 arises as the ratio of two appropriately
 * scaled chi-squared variates
 * X = (U1/d1) / (U2/d2)
 * U1, U2 ~ chi-squared distributions with d1 and d2 degrees of freedom
 * U1 and U2 are independent
 *
 * f(x;d1,d2)=1/B(d1/2,d2/2) * (d1/d2)^(d1/2) * x^(d1/2−1) * (1+d1/d2 x)^(−(d1+d2)/2) where x≥0
 */

int main() {
    double d1 = 10.0;
    double d2 = 12.0;

    /**
     * f(x;d1,d2)
     */
    cout << "density " << stats::df(1.5, d1, d2) << endl;  // 0.342627

    /**
     * P(X <= x)
     */
    cout << "cdf " << stats::pf(1.5, d1, d2) << endl;  // 0.75013

    /**
     * find x such that cdf(x) = a
     */
    cout << "quantile " << stats::qf(0.5, d1, d2) << endl;  // 0.98856
    cout << "quantile " << stats::qf(0.4, d1, d2) << endl;  // 0.844534
    cout << "quantile " << stats::qf(0.3, d1, d2) << endl;  // 0.712599

    /**
     * Random Sampling
     */
    cout << "single rand sample: " << stats::rf(d1, d2) << endl;
    std::vector<double> rand_samples = stats::rf<std::vector<double> >(10, 1, d1, d2);
    cout << "10 * 1 rand samples: " << rand_samples << endl;
}
