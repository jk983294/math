#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;

/**
 * recommended one is the BDCSVD class, which scale well for large problems and
 * automatically fall-back to the JacobiSVD class for smaller problems
 */

void basic_example();

int main() { basic_example(); }

/**
 * three observe point (1, 1), (1, 2), (2, 2)
 * fit with one line with least square ax + b = y
 * so resolve a, b as follow
 */
void basic_example() {
    MatrixXd A(3, 2);
    VectorXd b(3);
    A << 1, 1, 1, 1, 2, 1;  // A's second column is all 1, because b's coefficient is 1
    b << 1, 2, 2;
    cout << "Here is the matrix A:\n" << A << endl;
    cout << "Here is the right hand side b:\n" << b << endl;
    cout << "The least-squares solution is:\n" << A.bdcSvd(ComputeThinU | ComputeThinV).solve(b) << endl;
}
