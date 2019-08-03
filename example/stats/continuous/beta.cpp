#include <zerg_stl.h>
#include <iostream>
#include <stats.hpp>
#include <vector>

using namespace std;

/**
 * Beta distribution, continuous
 * f(x;α,β)=(1/B(α,β))*x^(α−1)*(1−x)^(β−1)
 */

int main() {
    double alpha = 3.0;
    double beta = 2.0;

    /**
     * f(x;α,β)
     */
    cout << "density " << stats::dbeta(0.5, alpha, beta) << endl;  // 1.5
    cout << "density " << stats::dbeta(0.5, alpha, beta) << endl;  // 1.5

    /**
     * P(X <= x)
     */
    cout << "cdf " << stats::pbeta(0.5, alpha, beta) << endl;  // 0.3125

    /**
     * find x such that cdf(x) = a
     */
    cout << "quantile " << stats::qbeta(0.5, alpha, beta) << endl;  // 0.614272
    cout << "quantile " << stats::qbeta(0.4, alpha, beta) << endl;  // 0.5555
    cout << "quantile " << stats::qbeta(0.3, alpha, beta) << endl;  // 0.491595

    /**
     * Random Sampling
     */
    cout << "single rand sample: " << stats::rbeta(alpha, beta) << endl;
    std::vector<double> rand_samples = stats::rbeta<std::vector<double> >(10, 1, alpha, beta);
    cout << "10 * 1 rand samples: " << rand_samples << endl;
}
