#include <zerg_stl.h>
#include <iostream>
#include <stats.hpp>
#include <vector>

using namespace std;

/**
 * Log Normal distribution, continuous
 * https://en.wikipedia.org/wiki/Log-normal_distribution
 * Y = ln(X) ~ N(μ,σ), then X = exp(Y), has a log-normal distribution
 *
 * The potential returns of a stock can be graphed in a normal distribution
 * The prices of the stock however can be graphed in a log-normal distribution
 *
 * f(x;μ,σ)=(1/x) (1/sqrt(2π)σ) exp(−(lnx−μ)^2/2σ^2)
 */

int main() {
    double mu = 0.0;
    double sigma = 1.0;

    /**
     * f(x;μ,σ)
     */
    cout << "density " << stats::dlnorm(1.0, mu, sigma) << endl;  // 0.398942
    cout << "density " << stats::dlnorm(0, mu, sigma) << endl;    // 0

    /**
     * P(X <= x)
     */
    cout << "cdf " << stats::plnorm(1.0, mu, sigma) << endl;  // 0.5
    cout << "cdf " << stats::plnorm(0, mu, sigma) << endl;    // 0

    /**
     * find x such that cdf(x) = a
     */
    cout << "quantile " << stats::qlnorm(0.841345, mu, sigma) << endl;  // 2.71828
    cout << "quantile " << stats::qlnorm(0.5, mu, sigma) << endl;       // 1
    cout << "quantile " << stats::qlnorm(0.4, mu, sigma) << endl;       // 0.776198
    cout << "quantile " << stats::qlnorm(0.3, mu, sigma) << endl;       // 0.59191

    /**
     * Random Sampling
     */
    cout << "single rand sample: " << stats::rlnorm(mu, sigma) << endl;
    std::vector<double> rand_samples = stats::rlnorm<std::vector<double> >(10, 1, mu, sigma);
    cout << "10 * 1 rand samples: " << rand_samples << endl;
}
