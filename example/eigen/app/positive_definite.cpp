#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;

/**
 * a positive definite matrix has properties:
 * 1) all eigen value > 0
 * 2) all sub-determinant > 0
 * 3) all pivots > 0
 * 4) x^T * A * x > 0
 *
 * for example:
 * [x y]^T * [1 3; 3 10] * [x y] = x^2 + 6xy + 10y^2 = (x + 3y)^2 + y^2 >= 0
 * it has min value at (0, 0)
 */

bool is_positive_definite(const MatrixXd& mat);

int main() {
    MatrixXd m(2, 2);
    m << 1, 3, 3, 10;

    cout << is_positive_definite(m) << endl;  // true
    m(1, 1) = 9;
    cout << is_positive_definite(m) << endl;  // false
    m(1, 1) = 8;
    cout << is_positive_definite(m) << endl;  // false
}

bool is_positive_definite(const MatrixXd& mat) {
    Eigen::LLT<Eigen::MatrixXd> lltOfA(mat);  // compute the Cholesky decomposition
    return lltOfA.info() != Eigen::NumericalIssue;
}
