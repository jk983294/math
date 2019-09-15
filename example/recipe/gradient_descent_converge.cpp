#include <zerg_stl.h>
#include <iostream>

using namespace std;

/**
 * minimize function using gradient descent
 *
 * f(x0, x1)= 0.5 (x0^2 + b * x1^2),  0 < b <= 1
 * gradient vector is <x0, bx1>
 *
 * eigen value: 1, b, when b -> 0, then it is hard to converge
 */

const double precision = 1e-9;
const double const_learning_rate = 0.1;

double f(double x0, double x1, double b) { return 0.5 * (x0 * x0 + b * x1 * x1); }

void update_ordinary_gradient(double &x0, double &x1, double b, double learning_rate) {
    x0 -= x0 * learning_rate;
    x1 -= b * x1 * learning_rate;
}

void test_ordinary_gradient(double b) {
    double x0 = b, x1 = 1;  // start point
    double f_current = f(x0, x1, b);
    int iter_count = 0;

    while (true) {
        update_ordinary_gradient(x0, x1, b, const_learning_rate);
        double f_new = f(x0, x1, b);
        double delta = f_new - f_current;
        if (abs(delta) <= precision) break;
        f_current = f_new;
        ++iter_count;
    }

    cout << "param b=" << b << ", after " << iter_count << " round, converge to <" << x0 << ", " << x1
         << "> value: " << f_current << endl;
}

void test_momentum_gradient(double b) {
    double x0 = b, x1 = 1;  // start point
    double f_current = f(x0, x1, b);
    double learning_rate = std::pow(2.0 / (1 + std::sqrt(b)), 2.0);
    double beta = std::pow((1 - std::sqrt(b)) / (1 + std::sqrt(b)), 2.0);
    double z0 = 0, z1 = 0;
    int iter_count = 0;

    while (true) {
        // z(k) = gradient(k) + beta * z(k - 1)
        z0 = x0 + beta * z0;
        z1 = b * x1 + beta * z1;
        // x(k+1) = x(k) - s * z(k)
        x0 -= learning_rate * z0;
        x1 -= learning_rate * z1;

        double f_new = f(x0, x1, b);
        double delta = f_new - f_current;
        if (abs(delta) <= precision) break;
        f_current = f_new;
        ++iter_count;
    }

    cout << "momentum gradient param b=" << b << ", after " << iter_count << " round, converge to <" << x0 << ", " << x1
         << "> value: " << f_current << endl;
}

int main() {
    test_ordinary_gradient(1);
    test_ordinary_gradient(0.1);
    test_ordinary_gradient(0.01);
    test_ordinary_gradient(0.001);

    test_momentum_gradient(1);
    test_momentum_gradient(0.1);
    test_momentum_gradient(0.01);
    test_momentum_gradient(0.001);
}
