#include <zerg_stl.h>
#include <iostream>
#include <stats.hpp>
#include <vector>

using namespace std;

/**
 * Uniform distribution, continuous
 *
 * f(x;a,b)=1/(b−a) where a≤x≤b
 */

int main() {
    double a = -1.0;
    double b = 1.0;

    /**
     * f(x;a,b)
     */
    cout << "density " << stats::dunif(1.0, a, b) << endl;  // 0.5
    cout << "density " << stats::dunif(0, a, b) << endl;    // 0.5

    /**
     * P(X <= x)
     */
    cout << "cdf " << stats::punif(0.5, a, b) << endl;  // 0.75
    cout << "cdf " << stats::punif(0, a, b) << endl;    // 0.5

    /**
     * find x such that cdf(x) = a
     */
    cout << "quantile " << stats::qunif(0.841345, a, b) << endl;  // 0.68269
    cout << "quantile " << stats::qunif(0.5, a, b) << endl;       // 0
    cout << "quantile " << stats::qunif(0.4, a, b) << endl;       // -0.2
    cout << "quantile " << stats::qunif(0.3, a, b) << endl;       // -0.4

    /**
     * Random Sampling
     */
    cout << "single rand sample: " << stats::runif(a, b) << endl;
    std::vector<double> rand_samples = stats::runif<std::vector<double> >(10, 1, a, b);
    cout << "10 * 1 rand samples: " << rand_samples << endl;
}
