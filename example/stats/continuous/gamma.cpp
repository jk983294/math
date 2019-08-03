#include <zerg_stl.h>
#include <iostream>
#include <stats.hpp>
#include <vector>

using namespace std;

/**
 * Gamma distribution, continuous
 * https://en.wikipedia.org/wiki/Gamma_distribution
 *
 * f(x;k,θ)= ( x^(k−1) * exp(−x/θ) ) / (θ^k * Γ(k)) where x>0
 */

int main() {
    double shape = 5.0;
    double scale = 4.0;

    /**
     * f(x;shape,scale)
     */
    cout << "density " << stats::dgamma(0.5, shape, scale) << endl;  // 2.24431e-06

    /**
     * P(X <= x)
     */
    cout << "cdf " << stats::pgamma(20.5, shape, scale) << endl;  // 0.581159

    /**
     * find x such that cdf(x) = a
     */
    cout << "quantile " << stats::qgamma(0.5, shape, scale) << endl;  // 18.6836
    cout << "quantile " << stats::qgamma(0.4, shape, scale) << endl;  // 16.5909
    cout << "quantile " << stats::qgamma(0.3, shape, scale) << endl;  // 14.5344

    /**
     * Random Sampling
     */
    cout << "single rand sample: " << stats::rgamma(shape, scale) << endl;
    std::vector<double> rand_samples = stats::rgamma<std::vector<double> >(10, 1, shape, scale);
    cout << "10 * 1 rand samples: " << rand_samples << endl;
}
