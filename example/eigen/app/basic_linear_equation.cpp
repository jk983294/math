#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;

/**
 * Ax=b
 */

void general_case();
void positive_definite();
void check_solution_exist();
void reuse();

int main() {
    general_case();
    positive_definite();
    check_solution_exist();
    reuse();
}

void general_case() {
    Matrix3f A;
    Vector3f b;
    A << 1, 2, 3, 4, 5, 6, 7, 8, 10;
    b << 3, 3, 4;
    cout << "Here is the matrix A:\n" << A << endl;
    cout << "Here is the vector b:\n" << b << endl;
    Vector3f x = A.colPivHouseholderQr().solve(b);
    cout << "The solution is:\n" << x << endl;
}

void positive_definite() {
    /**
     * for positive definite,  good choice is the LLT or LDLT decomposition
     * below also demo the b can be extended not only a vector but a matrix
     * so that each solution correspond to a column of b
     */
    Matrix2f A, b;
    A << 2, -1, -1, 3;
    b << 1, 2, 3, 1;
    cout << "Here is the matrix A:\n" << A << endl;
    cout << "Here is the right hand side b:\n" << b << endl;
    Matrix2f x = A.ldlt().solve(b);
    cout << "The solution is:\n" << x << endl;
}

/**
 * Only you know what error margin you want to allow for a solution to be considered valid
 */
void check_solution_exist() {
    MatrixXd A = MatrixXd::Random(100, 100);
    MatrixXd b = MatrixXd::Random(100, 1);
    MatrixXd x = A.fullPivLu().solve(b);
    double relative_error = (A * x - b).norm() / b.norm();  // norm() is L2 norm
    cout << "The relative error is:\n" << relative_error << endl;
}

void reuse() {
    Matrix2f A, b;
    LLT<Matrix2f> llt;
    A << 2, -1, -1, 3;
    b << 1, 2, 3, 1;
    cout << "Here is the matrix A:\n" << A << endl;
    cout << "Here is the right hand side b:\n" << b << endl;
    cout << "Computing LLT decomposition..." << endl;
    llt.compute(A);
    cout << "The solution is:\n" << llt.solve(b) << endl;
    A(1, 1)++;
    cout << "The matrix A is now:\n" << A << endl;
    cout << "Computing LLT decomposition..." << endl;
    llt.compute(A);
    cout << "The solution is now:\n" << llt.solve(b) << endl;
}
