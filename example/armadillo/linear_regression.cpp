#include <armadillo>
#include <iostream>
#include <stats.hpp>

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

    /**
     * solve不会改变原 A/b 的值
     */
    vec x = solve(A, b);  // throw if no solution found
    cout << A << endl;
    cout << b << endl;
    cout << x << endl;  // [2, 1, -1]^T

    vec coef;
    bool status = solve(coef, A, b);  // no throw
    cout << "status=" << status << endl;
    cout << A << endl;
    cout << b << endl;
    cout << "coef:\n" << coef << endl;  // [2, 1, -1]^T

    const arma::colvec& fitted = A * coef;
    const arma::colvec& residual = b - fitted;

    const int n = A.n_rows;
    const int k = A.n_cols;
    double s2 = arma::dot(residual, residual) / (n - k);
    const arma::colvec& std_err = arma::sqrt(s2 * arma::diagvec(arma::pinv(arma::trans(A) * A)));
    const arma::colvec& tvalue = coef / std_err;

    const arma::colvec& ytot = b - arma::mean(b);
    double r_squared = 1 - arma::dot(residual, residual) / arma::dot(ytot, ytot);

    cout << "fitted.values:\n" << fitted << endl;
    cout << "residuals:\n" << residual << endl;
    cout << "stderr:\n" << std_err << endl;
    cout << "tvalue:\n" << tvalue << endl;
    cout << "rsquared: " << r_squared << endl;

    for (size_t i = 0; i < tvalue.size(); ++i) {
        double p_value = 2 * (1 - stats::pt(std::abs(tvalue[i]), b.size() - tvalue.size()));
        printf("i=%zu, p_value=%f\n", i, p_value);
    }
    return 0;
}
