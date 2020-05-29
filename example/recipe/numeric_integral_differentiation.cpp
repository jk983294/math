#include <iostream>
#include <cmath>

using namespace std;

double func(double x) { return x * x * x; }
double func_integral(double a, double b) { return (std::pow(b, 4) - std::pow(a, 4)) / 4.0; }
double func_derative(double x) { return 3.0 * x * x; }
double func_second_order_derative(double x) { return 6.0 * x; }

double integral_midpoint(double a, double b, double delta) {
    long count = (b - a) / delta;
    double ret = 0;
    for (int i = 0; i < count; ++i) {
        double left = a + i * delta;
        double right = a + (i + 1) * delta;
        ret += func((left + right) / 2);
    }
    return ret * delta;
}

double integral_trapezoidal(double a, double b, double delta) {
    long count = (b - a) / delta;
    double ret = 0;
    for (int i = 0; i < count; ++i) {
        double left = a + i * delta;
        double right = a + (i + 1) * delta;
        ret += (func(left) + func(right)) / 2;
    }
    return ret * delta;
}

double integral_simpson(double a, double b, double delta) {
    long count = (b - a) / delta;
    double ret = 0;
    for (int i = 0; i < count; ++i) {
        double left = a + i * delta;
        double middle = a + (i + 1) * delta;
        double right = a + (i + 2) * delta;
        ret += (func(left) + 4 * func(middle) + func(right)) / 6;
    }
    return ret * delta;
}

double forward_difference(double x, double delta) {
    return (func(x + delta) - func(x)) / delta;
}

double backward_difference(double x, double delta) {
    return (func(x) - func(x - delta)) / delta;
}

double centered_difference(double x, double delta) {
    return (func(x + delta) - func(x - delta)) / (delta * 2);
}

double centered_second_order_difference(double x, double delta) {
    return (func(x + delta) - 2 * func(x) + func(x - delta)) / (delta * delta);
}

int main()
{
    double a = 0, b = 1.0;
    cout << "numerical integral" << endl;
    cout << "expected integral=" << func_integral(a, b) << endl;
    cout << "integral_midpoint delta=" << 1e-1 << " ret=" << integral_midpoint(a, b, 1e-1) << endl;
    cout << "integral_midpoint delta=" << 1e-2 << " ret=" << integral_midpoint(a, b, 1e-2) << endl;
    cout << "integral_midpoint delta=" << 1e-3 << " ret=" << integral_midpoint(a, b, 1e-3) << endl;
    cout << "integral_midpoint delta=" << 1e-4 << " ret=" << integral_midpoint(a, b, 1e-4) << endl;

    cout << "integral_trapezoidal delta=" << 1e-1 << " ret=" << integral_trapezoidal(a, b, 1e-1) << endl;
    cout << "integral_trapezoidal delta=" << 1e-2 << " ret=" << integral_trapezoidal(a, b, 1e-2) << endl;
    cout << "integral_trapezoidal delta=" << 1e-3 << " ret=" << integral_trapezoidal(a, b, 1e-3) << endl;
    cout << "integral_trapezoidal delta=" << 1e-4 << " ret=" << integral_trapezoidal(a, b, 1e-4) << endl;

    cout << "integral_simpson delta=" << 1e-1 << " ret=" << integral_simpson(a, b, 1e-1) << endl;
    cout << "integral_simpson delta=" << 1e-2 << " ret=" << integral_simpson(a, b, 1e-2) << endl;
    cout << "integral_simpson delta=" << 1e-3 << " ret=" << integral_simpson(a, b, 1e-3) << endl;
    cout << "integral_simpson delta=" << 1e-4 << " ret=" << integral_simpson(a, b, 1e-4) << endl;

    cout << "\nnumerical derative" << endl;
    cout << "expected derative=" << func_derative(b) << endl;
    cout << "forward_difference delta=" << 1e-1 << " ret=" << forward_difference(b, 1e-1) << endl;
    cout << "forward_difference delta=" << 1e-2 << " ret=" << forward_difference(b, 1e-2) << endl;
    cout << "forward_difference delta=" << 1e-3 << " ret=" << forward_difference(b, 1e-3) << endl;
    cout << "forward_difference delta=" << 1e-4 << " ret=" << forward_difference(b, 1e-4) << endl;
    cout << "backward_difference delta=" << 1e-1 << " ret=" << backward_difference(b, 1e-1) << endl;
    cout << "backward_difference delta=" << 1e-2 << " ret=" << backward_difference(b, 1e-2) << endl;
    cout << "backward_difference delta=" << 1e-3 << " ret=" << backward_difference(b, 1e-3) << endl;
    cout << "backward_difference delta=" << 1e-4 << " ret=" << backward_difference(b, 1e-4) << endl;
    cout << "centered_difference delta=" << 1e-1 << " ret=" << centered_difference(b, 1e-1) << endl;
    cout << "centered_difference delta=" << 1e-2 << " ret=" << centered_difference(b, 1e-2) << endl;
    cout << "centered_difference delta=" << 1e-3 << " ret=" << centered_difference(b, 1e-3) << endl;
    cout << "centered_difference delta=" << 1e-4 << " ret=" << centered_difference(b, 1e-4) << endl;

    cout << "\nexpected derative=" << func_second_order_derative(b) << endl;
    cout << "centered_second_order_difference delta=" << 1e-1 << " ret=" << centered_second_order_difference(b, 1e-1) << endl;
    cout << "centered_second_order_difference delta=" << 1e-2 << " ret=" << centered_second_order_difference(b, 1e-2) << endl;
    cout << "centered_second_order_difference delta=" << 1e-3 << " ret=" << centered_second_order_difference(b, 1e-3) << endl;
    cout << "centered_second_order_difference delta=" << 1e-4 << " ret=" << centered_second_order_difference(b, 1e-4) << endl;
    return 0;
}
