#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;

/**
 * for large matrix, inverse computations are often replaced by solve() operations,
 * and the determinant is often not a good way of checking if a matrix is invertible
 * for very small matrices, the above is not true, and inverse and determinant can be very useful.
 *
 * If matrix is of a very small fixed size (at most 4x4) this allows Eigen to avoid performing a LU decomposition,
 * and instead use formulas that are more efficient on such small matrices.
 */

int main() {
    Matrix3f A;
    A << 1, 2, 1, 2, 1, 0, -1, 1, 2;
    cout << "Here is the matrix A:\n" << A << endl;
    cout << "The determinant of A is " << A.determinant() << endl;
    cout << "The inverse of A is:\n" << A.inverse() << endl;
}
