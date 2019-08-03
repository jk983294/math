#include <zerg_stl.h>
#include <iostream>
#include <stats.hpp>
#include <vector>

using namespace std;

/**
 * Exponential distribution, continuous
 * the time between events in a Poisson point process
 * https://en.wikipedia.org/wiki/Exponential_distribution
 *
 * f(x;λ)= λ exp(−λx) where x≥0
 */

int main() {
    double rate = 2.0;  // average 2 unit time, one event occur

    /**
     * f(x;λ)
     */
    cout << "density " << stats::dexp(0.5, rate) << endl;  // 0.735759

    /**
     * P(X <= x)
     */
    cout << "cdf " << stats::pexp(0.5, rate) << endl;  // 0.632121

    /**
     * find x such that cdf(x) = a
     */
    cout << "quantile " << stats::qexp(0.5, rate) << endl;  // 0.346574
    cout << "quantile " << stats::qexp(0.4, rate) << endl;  // 0.255413
    cout << "quantile " << stats::qexp(0.3, rate) << endl;  // 0.178337

    /**
     * Random Sampling
     */
    cout << "single rand sample: " << stats::rexp(rate) << endl;
    std::vector<double> rand_samples = stats::rexp<std::vector<double> >(10, 1, rate);
    cout << "10 * 1 rand samples: " << rand_samples << endl;
}
