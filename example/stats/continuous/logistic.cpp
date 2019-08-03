#include <zerg_stl.h>
#include <iostream>
#include <stats.hpp>
#include <vector>

using namespace std;

/**
 * Logistic distribution, continuous
 * https://en.wikipedia.org/wiki/Logistic_distribution
 *
 * It resembles the normal distribution in shape but has heavier tails (higher kurtosis)
 *
 * f(x;μ,σ)=exp(−(x−μ)/σ) / (σ*(1+exp(−(x−μ)/σ))^2)
 */

int main() {
    double mu = 0.0;
    double sigma = 1.0;

    /**
     * f(x;μ,σ)
     */
    cout << "density " << stats::dlogis(1.0, mu, sigma) << endl;  // 0.196612
    cout << "density " << stats::dlogis(0, mu, sigma) << endl;    // 0.25

    /**
     * P(X <= x)
     */
    cout << "cdf " << stats::plogis(1.0, mu, sigma) << endl;  // 0.731059
    cout << "cdf " << stats::plogis(0, mu, sigma) << endl;    // 0.5

    /**
     * find x such that cdf(x) = a
     */
    cout << "quantile " << stats::qlogis(0.841345, mu, sigma) << endl;  // 1.66827
    cout << "quantile " << stats::qlogis(0.5, mu, sigma) << endl;       // 0
    cout << "quantile " << stats::qlogis(0.4, mu, sigma) << endl;       // -0.405465
    cout << "quantile " << stats::qlogis(0.3, mu, sigma) << endl;       // -0.847298

    /**
     * Random Sampling
     */
    cout << "single rand sample: " << stats::rlogis(mu, sigma) << endl;
    std::vector<double> rand_samples = stats::rlogis<std::vector<double> >(10, 1, mu, sigma);
    cout << "10 * 1 rand samples: " << rand_samples << endl;
}
