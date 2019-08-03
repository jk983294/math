#include <zerg_stl.h>
#include <iostream>
#include <stats.hpp>
#include <vector>

using namespace std;

/**
 * Cauchy distribution, continuous
 * https://en.wikipedia.org/wiki/Cauchy_distribution
 * f(x;μ,σ)= 1 / ( πσ[1+(x−μσ)^2] )
 */

int main() {
    double mu = 0.0;
    double sigma = 1.0;

    /**
     * f(x;α,β)
     */
    cout << "density " << stats::dcauchy(0.5, mu, sigma) << endl;  // 0.254648
    cout << "density " << stats::dcauchy(0.8, mu, sigma) << endl;  // 0.194091

    /**
     * P(X <= x)
     */
    cout << "cdf " << stats::pcauchy(0.5, mu, sigma) << endl;  // 0.647584

    /**
     * find x such that cdf(x) = a
     */
    cout << "quantile " << stats::qcauchy(0.647584, mu, sigma) << endl;  // 0.500002
    cout << "quantile " << stats::qcauchy(0.5, mu, sigma) << endl;       // 0
    cout << "quantile " << stats::qcauchy(0.4, mu, sigma) << endl;       // -0.32492
    cout << "quantile " << stats::qcauchy(0.3, mu, sigma) << endl;       // -0.726543

    /**
     * Random Sampling
     */
    cout << "single rand sample: " << stats::rcauchy(mu, sigma) << endl;
    std::vector<double> rand_samples = stats::rcauchy<std::vector<double> >(10, 1, mu, sigma);
    cout << "10 * 1 rand samples: " << rand_samples << endl;
}
