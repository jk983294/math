#include <zerg_stl.h>
#include <iostream>
#include <stats.hpp>
#include <vector>

using namespace std;

/**
 * Laplace distribution, continuous
 * https://en.wikipedia.org/wiki/Laplace_distribution
 * The difference between two independent identically distributed exponential random variables
 * is governed by a Laplace distribution
 *
 * f(x;μ,σ)=(1/2σ)exp(−|x−μ|/σ)
 */

int main() {
    double mu = 0.0;
    double sigma = 1.0;

    /**
     * f(x;μ,σ)
     */
    cout << "density " << stats::dlaplace(0.5, mu, sigma) << endl;  // 0.303265
    cout << "density " << stats::dlaplace(0.8, mu, sigma) << endl;  // 0.224664

    /**
     * P(X <= x)
     */
    cout << "cdf " << stats::plaplace(0.5, mu, sigma) << endl;  // 0.696735

    /**
     * find x such that cdf(x) = a
     */
    cout << "quantile " << stats::qlaplace(0.647584, mu, sigma) << endl;  // 0.349796
    cout << "quantile " << stats::qlaplace(0.5, mu, sigma) << endl;       // 0
    cout << "quantile " << stats::qlaplace(0.4, mu, sigma) << endl;       // -0.223144
    cout << "quantile " << stats::qlaplace(0.3, mu, sigma) << endl;       // -0.510826

    /**
     * Random Sampling
     */
    cout << "single rand sample: " << stats::rlaplace(mu, sigma) << endl;
    std::vector<double> rand_samples = stats::rlaplace<std::vector<double> >(10, 1, mu, sigma);
    cout << "10 * 1 rand samples: " << rand_samples << endl;
}
