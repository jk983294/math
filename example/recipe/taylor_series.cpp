#include <cmath>
#include <iostream>

using namespace std;

/**
 * approximation of ln(1 + x)
 */
double taylor_series(double x, int n) {
    double ret = 0;
    for (int i = 1; i < n; i++) {
        double val = std::pow(x, i) / i;
        if (i % 2 == 0) val = -val;
        ret += val;
    }
    return ret;
}

int main() {
    double x = 1.5;  // not in convergent area
    std::cout << taylor_series(x, 20) << " <--> " << std::log(1 + x) << std::endl;
    x = 0.5;  // in convergent area
    std::cout << taylor_series(x, 20) << " <--> " << std::log(1 + x) << std::endl;
    return 0;
}
