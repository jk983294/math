#include <zerg_stl.h>
#include <iostream>
#include <stats.hpp>
#include <vector>

using namespace std;

/**
 * Gaussian distribution, continuous
 * https://en.wikipedia.org/wiki/Normal_distribution
 *
 * f(x;μ,σ)=(1/sqrt(2π)σ) exp(−(x−μ)^2/2σ^2)
 */

int main() {
    double mu = 0.0;
    double sigma = 1.0;

    /**
     * f(x;μ,σ)
     */
    cout << "density " << stats::dnorm(1.0, mu, sigma) << endl;  // 0.241971
    cout << "density " << stats::dnorm(0, mu, sigma) << endl;    // 0.398942

    /**
     * P(X <= x)
     */
    cout << "cdf " << stats::pnorm(1.0, mu, sigma) << endl;  // 0.841345
    cout << "cdf " << stats::pnorm(0, mu, sigma) << endl;    // 0.5

    /**
     * find x such that cdf(x) = a
     */
    cout << "quantile " << stats::qnorm(0.841345, mu, sigma) << endl;  // 1
    cout << "quantile " << stats::qnorm(0.5, mu, sigma) << endl;       // 0
    cout << "quantile " << stats::qnorm(0.4, mu, sigma) << endl;       // -0.253347
    cout << "quantile " << stats::qnorm(0.3, mu, sigma) << endl;       // -0.524401

    /**
     * Random Sampling
     */
    cout << "single rand sample: " << stats::rnorm(mu, sigma) << endl;
    std::vector<double> rand_samples = stats::rnorm<std::vector<double> >(10, 1, mu, sigma);
    cout << "10 * 1 rand samples: " << rand_samples << endl;
}
