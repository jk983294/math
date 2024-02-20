#include <zerg_stl.h>
#include <iostream>
#include <stats.hpp>
#include <vector>
#include <math_dist.h>

using namespace std;
using namespace ornate;

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
    printf("density %f,%f\n", stats::dnorm(1.0, mu, sigma), norm_pdf(1.0, mu, sigma)); // 0.241971
    printf("density %f,%f\n", stats::dnorm(0, mu, sigma), norm_pdf(0, mu, sigma)); // 0.398942
    printf("density %f,%f\n", stats::dnorm(-1., mu, sigma), norm_pdf(-1., mu, sigma)); // 0.241971

    /**
     * P(X <= x)
     */
    printf("cdf %f,%f\n", stats::pnorm(1.0, mu, sigma), norm_cdf(1.0, mu, sigma)); // 0.841345
    printf("cdf %f,%f\n", stats::pnorm(0, mu, sigma), norm_cdf(0, mu, sigma)); // 0.5
    printf("cdf %f,%f\n", stats::pnorm(-1., mu, sigma), norm_cdf(-1., mu, sigma)); // 0.158655

    /**
     * find x such that cdf(x) = a
     */
    printf("icdf %f,%f\n", stats::qnorm(0.841345, mu, sigma), norm_icdf(0.841345, mu, sigma)); // 1
    printf("icdf %f,%f\n", stats::qnorm(0.5, mu, sigma), norm_icdf(0.5, mu, sigma)); // 0
    printf("icdf %f,%f\n", stats::qnorm(0.4, mu, sigma), norm_icdf(0.4, mu, sigma)); // -0.253347
    printf("icdf %f,%f\n", stats::qnorm(0.3, mu, sigma), norm_icdf(0.3, mu, sigma)); // -0.524401

    /**
     * Random Sampling
     */
    cout << "single rand sample: " << stats::rnorm(mu, sigma) << endl;
    std::vector<double> rand_samples = stats::rnorm<std::vector<double> >(10, 1, mu, sigma);
    cout << "10 * 1 rand samples: " << rand_samples << endl;
}
