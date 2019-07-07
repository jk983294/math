#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;

/**
 * Certain decompositions are able to compute the rank of a matrix
 */

void basic_example();
void change_threshold();

int main() {
    basic_example();
    change_threshold();
}

void basic_example() {
    Matrix3f A;
    A << 1, 2, 5, 2, 1, 4, 3, 0, 3;
    cout << "Here is the matrix A:\n" << A << endl;
    FullPivLU<Matrix3f> lu_decomp(A);
    cout << "The rank of A is " << lu_decomp.rank() << endl;
    cout << "Here is a matrix whose columns form a basis of the null-space of A:\n" << lu_decomp.kernel() << endl;
    cout << "Here is a matrix whose columns form a basis of the column-space of A:\n"
         << lu_decomp.image(A) << endl;  // yes, have to pass the original A
}

/**
 * any rank computation depends on the choice of threshold,
 * since practically no floating-point matrix is exactly rank-deficient
 *
 * only you know what is the right threshold for your application.
 *
 * The decomposition itself, i.e. the compute() method, is independent of the threshold.
 * You don't need to recompute the decomposition after you've changed the threshold.
 */
void change_threshold() {
    Matrix2d A;
    A << 2, 1, 2, 0.9999999999;
    FullPivLU<Matrix2d> lu(A);
    cout << "By default, the rank of A is found to be " << lu.rank() << endl;
    lu.setThreshold(1e-5);
    cout << "With threshold 1e-5, the rank of A is found to be " << lu.rank() << endl;
}
