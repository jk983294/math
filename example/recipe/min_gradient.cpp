#include <zerg_stl.h>
#include <iostream>

using namespace std;

/**
 * minimize function using gradient descent
 *
 * f(x, y)= x^2 + y^2
 * gradient vector is <2x, 2y>, its direction is point to greater f(x, y)
 * so we choose <-2x, -2y> to update new position
 *
 * as close to min position, the update rate will be smaller as derivative get more and more flat
 */

const double precision = 1e-9;
/**
 * carefully choose learning rate, if you choose 1.0, it will not converge
 */
const double learning_rate = 0.3;

double f(double x, double y) { return x * x + y * y; }

void update_along_gradient(double& x, double& y) {
    x += -2.0 * x * learning_rate;
    y += -2.0 * y * learning_rate;
}

int main() {
    double x_pos = 10, y_pos = 10;  // start point
    double f_current = f(x_pos, y_pos);
    int iter_count = 0;

    while (true) {
        update_along_gradient(x_pos, y_pos);
        double f_new = f(x_pos, y_pos);
        double delta = f_new - f_current;
        if (abs(delta) <= precision) break;
        f_current = f_new;
        ++iter_count;
        cout << x_pos << " " << y_pos << " " << delta << endl;
    }

    cout << "after " << iter_count << " round, converge to <" << x_pos << ", " << y_pos << "> value: " << f_current
         << endl;
}
