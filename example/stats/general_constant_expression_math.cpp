#include <gcem.hpp>
#include <iostream>

using namespace std;

/**
 * compile time calculation
 */

void algorithm();
void basic();
void special();

int main() {
    algorithm();
    basic();
    special();
}

void special() {
    cout << "binomial coefficient choose C(3, 2) = " << gcem::binomial_coef(3, 2) << endl;  // 3
    cout << "binomial coefficient choose C(4, 2) = " << gcem::binomial_coef(4, 2) << endl;  // 6

    cout << "beta " << gcem::beta(2, 3) << endl;  // 0.0833333
    cout << "beta " << gcem::beta(2, 4) << endl;  // 0.05
}

void basic() {
    cout << "exp(x)−1 expm1(1) = " << gcem::expm1(1) << endl;          // 1.71828
    cout << "exp(x)−1 expm1(2) = " << gcem::expm1(2) << endl;          // 6.38906
    cout << "loge (x+1) log1p(1) = " << gcem::log1p(1.71828) << endl;  // 1
    cout << "loge (x+1) log1p(2) = " << gcem::log1p(6.38906) << endl;  // 2

    cout << "factorial int factorial(3) = " << gcem::factorial(3) << endl;                     // 6
    cout << "factorial double factorial(x) = tgamma(x+1) = " << gcem::factorial(3.1) << endl;  // 6.81262

    cout << "sign function sgn(-4.2) = " << gcem::sgn(-4.2) << endl;  // -1
    cout << "sign function sgn(0) = " << gcem::sgn(0) << endl;        // 0
    cout << "sign function sgn(42) = " << gcem::sgn(42) << endl;      // 1

    cout << "trunc(-4.6) = " << gcem::trunc(-4.6) << endl;  // -4
    cout << "trunc(1) = " << gcem::trunc(1) << endl;        // 1
    cout << "trunc(4.6) = " << gcem::trunc(4.6) << endl;    // 4
}

void algorithm() {
    cout << "greatest common divisor (12, 18) = " << gcem::gcd(12, 18) << endl;  // 6
    cout << "greatest common divisor (12, 3)  = " << gcem::gcd(12, 3) << endl;   // 3
    cout << "greatest common divisor (7, 2)   = " << gcem::gcd(7, 2) << endl;    // 1

    cout << "least common multiple (12, 18) = " << gcem::lcm(12, 18) << endl;  // 36
    cout << "least common multiple (12, 3)  = " << gcem::lcm(12, 3) << endl;   // 12
    cout << "least common multiple  (7, 2)  = " << gcem::lcm(7, 2) << endl;    // 14
}
