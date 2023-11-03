#include <Eigen/Dense>
#include <iostream>
#include <stats.hpp>

using namespace std;
using namespace Eigen;

/**
 * Y = A(0)*X[0] + ... + A(N-1)*X[N-1] + A(N)
 * train data generated from y= 2*x[0] + x[1] -1
 */

int main() {
    MatrixXd A;
    A.resize(9, 3);
    // col major
    std::vector<double> a_data = {0, 1, 2, 0, 0, 1, 1, 2, 2, 0, 0, 0, 1, 2, 1, 2, 1., 2, 1, 1., 1, 1, 1., 1, 1, 1, 1};
    // A <<
    // 0, 0, 1,
    // 1, 0, 1,
    // 2, 0, 1,
    // 0, 1, 1,
    // 0, 2, 1,
    // 1, 1., 1,
    // 1, 2., 1,
    // 2, 1., 1,
    // 2, 2., 1;  // << operator use row major
    double* pA = A.data();
    std::copy(a_data.begin(), a_data.end(), pA);
    VectorXd b;
    b.resize(9);
    std::vector<double> b_data = {-1.0, 1.0, 3.0, 0.0, 1.0, 2.0, 3.0, 4, 5};
    double* pb= b.data();
    std::copy(b_data.begin(), b_data.end(), pb);

    VectorXd coef = A.bdcSvd(ComputeThinU | ComputeThinV).solve(b);
    cout << "A:\n" << A << endl;
    cout << "b:\n" << b << endl;
    cout << "coef:\n" << coef << endl;  // [2, 1, -1]^T

    const VectorXd& fitted = A * coef;
    const VectorXd& residual = b - fitted;

    const int n = A.rows();
    const int k = A.cols();
    double s2 = residual.dot(residual) / (n - k);
    auto& pinv = (A.transpose() * A).completeOrthogonalDecomposition().pseudoInverse();
    const VectorXd& std_err = (s2 * pinv.diagonal()).array().sqrt();
    const VectorXd& tvalue = coef.array() / std_err.array();

    const VectorXd& ytot = b.array() - b.array().mean();
    double r_squared = 1 - residual.dot(residual) / ytot.dot(ytot);

    cout << "fitted.values:\n" << fitted << endl;
    cout << "residuals:\n" << residual << endl;
    cout << "stderr:\n" << std_err << endl;
    cout << "tvalue:\n" << tvalue << endl;
    cout << "rsquared: " << r_squared << endl;

    for (int i = 0; i < tvalue.size(); ++i) {
        double p_value = 2 * (1 - stats::pt(std::abs(tvalue[i]), b.size() - tvalue.size()));
        printf("i=%d, p_value=%f\n", i, p_value);
    }
    return 0;
}
