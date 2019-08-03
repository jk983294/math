#include <zerg_stl.h>
#include <iostream>
#include <stats.hpp>
#include <vector>

using namespace std;

/**
 * Chi-squared distribution, continuous
 * the distribution of a sum of the squares of k independent standard normal random variables
 * https://en.wikipedia.org/wiki/Chi-squared_distribution
 *
 * f(x;k)=x^(k/2−1)exp(−x/2) / (2^(k/2) Γ(k/2)) where x>0
 */

int main() {
    double k = 14.0;

    /**
     * f(x;k)
     */
    cout << "density " << stats::dchisq(4.5, k) << endl;  // 0.00949664

    /**
     * P(X <= x)
     */
    cout << "cdf " << stats::pchisq(4.5, k) << endl;  // 0.00837207

    /**
     * find x such that cdf(x) = a
     */
    cout << "quantile " << stats::qchisq(0.5, k) << endl;  // 13.3393
    cout << "quantile " << stats::qchisq(0.4, k) << endl;  // 12.0785
    cout << "quantile " << stats::qchisq(0.3, k) << endl;  // 10.8215

    /**
     * Random Sampling
     */
    cout << "single rand sample: " << stats::rchisq(k) << endl;
    std::vector<double> rand_samples = stats::rchisq<std::vector<double> >(10, 1, k);
    cout << "10 * 1 rand samples: " << rand_samples << endl;
}
