#include <zerg_stl.h>
#include <iostream>
#include <stats.hpp>
#include <vector>

using namespace std;

/**
 * Bernoulli distribution, discrete
 * f(x;p)= p^x * (1−p)^(1−x)
 */

int main() {
    double prob_param = 0.6;

    /**
     * P(X = x)
     */
    cout << "density " << stats::dbern(0, prob_param) << endl;  // 0.4
    cout << "density " << stats::dbern(1, prob_param) << endl;  // 0.6

    /**
     * P(X <= x)
     */
    cout << "cdf " << stats::pbern(0, prob_param) << endl;  // 0.4
    cout << "cdf " << stats::pbern(1, prob_param) << endl;  // 1.0

    /**
     * find x such that cdf(x) = a
     */
    cout << "quantile " << stats::qbern(0.5, prob_param) << endl;  // 1
    cout << "quantile " << stats::qbern(0.4, prob_param) << endl;  // 0
    cout << "quantile " << stats::qbern(0.3, prob_param) << endl;  // 0

    /**
     * Random Sampling
     */
    cout << "single rand sample: " << stats::rbern(prob_param) << endl;
    std::vector<int> rand_samples = stats::rbern<std::vector<int> >(10, 1, prob_param);
    cout << "10 * 1 rand samples: " << rand_samples << endl;
}
