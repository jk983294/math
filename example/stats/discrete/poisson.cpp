#include <zerg_stl.h>
#include <iostream>
#include <stats.hpp>
#include <vector>

using namespace std;

/**
 * Poisson distribution, continuous
 * the probability of a given number of events occurring in a fixed interval of time
 * if these events occur with a known constant rate and independently of the time since the last event
 *
 * https://en.wikipedia.org/wiki/Poisson_distribution
 *
 * f(x;λ)= λ^x * exp(−λ) / (x!) where x≥0
 */

int main() {
    double rate = 2.0;  // average 2 unit time, one event occur

    /**
     * f(x;λ)
     */
    cout << "density " << stats::dpois(2, rate) << endl;  // 0.270671

    /**
     * P(X <= x)
     */
    cout << "cdf " << stats::ppois(2, rate) << endl;  // 0.676676

    /**
     * find x such that cdf(x) = a
     */
    cout << "quantile " << stats::qpois(0.5, rate) << endl;  // 2
    cout << "quantile " << stats::qpois(0.4, rate) << endl;  // 1
    cout << "quantile " << stats::qpois(0.3, rate) << endl;  // 1

    /**
     * Random Sampling
     */
    cout << "single rand sample: " << stats::rpois(rate) << endl;
    std::vector<double> rand_samples = stats::rpois<std::vector<double> >(10, 1, rate);
    cout << "10 * 1 rand samples: " << rand_samples << endl;
}
