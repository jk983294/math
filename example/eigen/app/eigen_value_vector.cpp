#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;

int main() {
    Matrix2f A;
    A << 1, 2, 2, 3;
    cout << "Here is the matrix A:\n" << A << endl;
    SelfAdjointEigenSolver<Matrix2f> eigenSolver(A);
    if (eigenSolver.info() != Success) {
        cerr << "this matrix eigen solve does not converge" << endl;
        abort();
    }
    cout << "The eigenvalues of A are:\n" << eigenSolver.eigenvalues() << endl;
    cout << "Here's a matrix whose columns are eigen vectors of A \n"
         << "corresponding to these eigen values:\n"
         << eigenSolver.eigenvectors() << endl;
}
