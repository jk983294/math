#include <zerg_stl.h>
#include <iostream>
#include <stats.hpp>
#include <vector>

using namespace std;

/**
 * Binomial distribution, discrete
 * f(x;n,p)= C(n,x) p^x (1−p)^(n−x)
 */

int main() {
    double prob_param = 0.6;
    unsigned int n_trials = 4;

    /**
     * P(X = x)
     */
    cout << "density " << stats::dbinom(0, n_trials, prob_param) << endl;  // 0.0256
    cout << "density " << stats::dbinom(1, n_trials, prob_param) << endl;  // 0.1536

    /**
     * P(X <= x)
     */
    cout << "cdf " << stats::pbinom(0, n_trials, prob_param) << endl;  // 0.0256
    cout << "cdf " << stats::pbinom(1, n_trials, prob_param) << endl;  // 0.1792

    /**
     * find x such that cdf(x) = a
     */
    cout << "quantile " << stats::qbinom(0.5, n_trials, prob_param) << endl;     // 2
    cout << "quantile " << stats::qbinom(0.1536, n_trials, prob_param) << endl;  // 1
    cout << "quantile " << stats::qbinom(0.0255, n_trials, prob_param) << endl;  // 0

    /**
     * Random Sampling
     */
    cout << "single rand sample: " << stats::rbinom(n_trials, prob_param) << endl;
    std::vector<int> rand_samples = stats::rbinom<std::vector<int> >(10, 1, n_trials, prob_param);
    cout << "10 * 1 rand samples: " << rand_samples << endl;
}
