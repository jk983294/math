#include <zerg_stl.h>
#include <iostream>
#include <stats.hpp>
#include <vector>

using namespace std;

/**
 * Inverse Gamma distribution, continuous
 * https://en.wikipedia.org/wiki/Inverse-invgamma_distribution
 *
 * f(x;α,β)=(β^α/Γ(α)) * x^(−α−1) * exp(−β/x) where x>0
 */

int main() {
    double shape = 5.0;
    double scale = 4.0;

    /**
     * f(x;shape,scale)
     */
    cout << "density " << stats::dinvgamma(0.5, shape, scale) << endl;  // 0.916037

    /**
     * P(X <= x)
     */
    cout << "cdf " << stats::pinvgamma(0.8, shape, scale) << endl;  // 0.440493

    /**
     * find x such that cdf(x) = a
     */
    cout << "quantile " << stats::qinvgamma(0.5, shape, scale) << endl;  // 0.856364
    cout << "quantile " << stats::qinvgamma(0.4, shape, scale) << endl;  // 0.763852
    cout << "quantile " << stats::qinvgamma(0.3, shape, scale) << endl;  // 0.679075

    /**
     * Random Sampling
     */
    cout << "single rand sample: " << stats::rinvgamma(shape, scale) << endl;
    std::vector<double> rand_samples = stats::rinvgamma<std::vector<double> >(10, 1, shape, scale);
    cout << "10 * 1 rand samples: " << rand_samples << endl;
}
