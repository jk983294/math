#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;

/**
 * Y = A(0)*X[0] + ... + A(N-1)*X[N-1] + A(N)
 * train data generated from y= 2*x[0] + x[1] -1
 */

int main() {
    mat A;
    A.set_size(9, 2);
    A << 0 << 0 << 1 << endr << 1 << 0 << 1 << endr << 2 << 0 << 1 << endr << 0 << 1 << 1 << endr << 0 << 2 << 1 << endr
      << 1 << 1. << 1 << endr << 1 << 2. << 1 << endr << 2 << 1. << 1 << endr << 2 << 2. << 1 << endr;
    vec b;
    b.resize(9);
    b[0] = -1.0;
    b[1] = 1.0;
    b[2] = 3.0;
    b[3] = 0.0;
    b[4] = 1.0;
    b[5] = 2.0;
    b[6] = 3.0;
    b[7] = 4.0;
    b[8] = 5.0;

    vec x = solve(A, b);  // throw if no solution found
    cout << A << endl;
    cout << b << endl;
    cout << x << endl;

    vec x2;
    bool status = solve(x2, A, b);  // no throw
    cout << "status=" << status << endl;
    cout << A << endl;
    cout << b << endl;
    cout << x2 << endl;
    return 0;
}
