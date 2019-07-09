#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;

/**
 * a markov matrix has properties:
 * 1) all element >= 0
 * 2) sum of each column = 1
 * 3) eigen value is 1 which respond to steady state, the rest eigen values, |lamba(i)| < 1
 * 4) eigen value of A = eigen value of A^T
 *
 * y(t+1) = M * y(t)
 */

VectorXd naive_multiply_method(const MatrixXd& mat, const VectorXd& y, int times);
VectorXd eigen_value_method(const MatrixXd& mat, const VectorXd& y);
VectorXd eigen_value_1_method(const MatrixXd& mat, const VectorXd& y);

int main() {
    MatrixXd m(2, 2);
    m << 0.9, 0.2, 0.1, 0.8;
    VectorXd y(2);
    y << 0, 1000;  // init state

    VectorXd result100 = naive_multiply_method(m, y, 20);
    cout << "result:\n" << result100 << endl;

    VectorXd result = eigen_value_method(m, y);     // the accuracy is not good
    cout << "result:\n" << result << endl;

    VectorXd result1 = eigen_value_1_method(m, y);  // this solution is best
    cout << "result:\n" << result1 << endl;
}

VectorXd naive_multiply_method(const MatrixXd& mat, const VectorXd& y, int times) {
    VectorXd result = y;
    for (int i = 0; i < times; ++i) {
        VectorXd tmp = mat * result;
        result = tmp;
    }
    return result;
}

VectorXd eigen_value_method(const MatrixXd& mat, const VectorXd& y) {
    SelfAdjointEigenSolver<MatrixXd> eigenSolver(mat);
    if (eigenSolver.info() != Success) {
        cerr << "this markov matrix eigen solve does not converge" << endl;
        abort();
    }

    const VectorXd& eigens = eigenSolver.eigenvalues();
    MatrixXf::Index maxIndex;
    double maxEigen = eigens.maxCoeff(&maxIndex);
    cout << "The eigenvalues of A are:\n" << eigens << endl;
    cout << "the max eigen value of A are: " << maxEigen << " and its pos is " << maxIndex << endl;
    VectorXd distribution = eigenSolver.eigenvectors().col(maxIndex);
    cout << "corresponding to these eigen vector:\n" << distribution << endl;
    return (distribution.array() / distribution.array().sum()) * y.array().sum();
}

VectorXd eigen_value_1_method(const MatrixXd& mat, const VectorXd& y) {
    MatrixXd minus1mat = mat - MatrixXd::Identity(mat.rows(), mat.cols()); // |A - lambda * I|
    FullPivLU<MatrixXd> lu_decomp(minus1mat);
    VectorXd distribution = lu_decomp.kernel().col(0);
    cout << "eigen value 1's eigen vector:\n" << distribution << endl;
    return (distribution.array() / distribution.array().sum()) * y.array().sum();
}
