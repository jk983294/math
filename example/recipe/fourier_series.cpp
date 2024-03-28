#include <cmath>
#include <iostream>

using namespace std;

/**
 * y = x, x belong to [-pi, pi)
 */
double fourier_series(double x, int n) {
    double ret = 0;
    for (int i = 1; i < n; i++) {
        double val = 2. * std::sin(i * x) / i;
        if (i % 2 == 0) val = -val;
        ret += val;
    }
    return ret;
}

int main() {
    double x = 3.14 / 2;
    for (size_t i = 1; i < 100; i++) {
        double approx = fourier_series(x, i);
        double diff = std::abs(approx - x);
        std::cout << approx << " <--> " << x << " diff: " << diff << std::endl;
    }
    x = -3.14 / 2;
    std::cout << fourier_series(x, 100) << " <--> " << x << std::endl;
    return 0;
}
