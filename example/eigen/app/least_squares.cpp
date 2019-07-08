#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;

/**
 * An overdetermined system of equations, say Ax = b, has no solutions. In this case,
 * it makes sense to search for the vector x which is closest to being a solution, in the sense that
 * the difference Ax - b is as small as possible. This x is called the least square solution
 *
 * recommended one is the BDCSVD class, which scale well for large problems and
 * automatically fall-back to the JacobiSVD class for smaller problems
 *
 * the SVD decomposition is generally the most accurate but the slowest,
 * normal equations is the fastest but least accurate, and the QR decomposition is in between
 */

void basic_example();
void qr_example();
void normal_equation();

int main() {
    basic_example();
    qr_example();
    normal_equation();
}

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

void qr_example() {
    MatrixXd A(3, 2);
    VectorXd b(3);
    A << 1, 1, 1, 1, 2, 1;  // A's second column is all 1, because b's coefficient is 1
    b << 1, 2, 2;
    cout << "The solution using the QR decomposition is:\n" << A.colPivHouseholderQr().solve(b) << endl;
}

void normal_equation() {
    MatrixXd A(3, 2);
    VectorXd b(3);
    A << 1, 1, 1, 1, 2, 1;  // A's second column is all 1, because b's coefficient is 1
    b << 1, 2, 2;
    cout << "The solution using normal equations is:\n" << (A.transpose() * A).ldlt().solve(A.transpose() * b) << endl;
}
