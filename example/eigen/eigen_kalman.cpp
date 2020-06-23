#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;

class Kalman {
private:
    MatrixXd A;  // state transition matrix
    MatrixXd H;  // measurement matrix
    MatrixXd Q, R;
    MatrixXd xest, pest;  // state conditions
    double dt;

public:
    // constructor, sets up data structures
    Kalman() : dt(1.0) {
        A.setIdentity(6, 6);
        A(0, 2) = A(1, 3) = A(2, 4) = A(3, 5) = dt;

        H.setZero(2, 6);
        H(0, 0) = H(1, 1) = 1.0;

        Q.setIdentity(6, 6);
        R = 1000 * R.Identity(2, 2);

        xest.setZero(6, 1);
        pest.setZero(6, 6);
    }

    // sole member function: estimate model
    MatrixXd estimate(const MatrixXd& Z) {
        size_t n = Z.rows(), k = Z.cols();
        MatrixXd Y = MatrixXd::Zero(n, k);
        MatrixXd xprd, pprd, S, B, kalmangain;
        VectorXd y;

        for (unsigned int i = 0; i < n; i++) {
            VectorXd z = Z.row(i).transpose();
            // predicted state and covariance
            xprd = A * xest;
            pprd = A * pest * A.transpose() + Q;

            // estimation
            S = H * pprd.transpose() * H.transpose() + R;
            B = H * pprd.transpose();

            kalmangain = S.ldlt().solve(B).transpose();

            // estimated state and covariance
            xest = xprd + kalmangain * (z - H * xprd);
            pest = pprd - kalmangain * H * pprd;

            // compute the estimated measurements
            y = H * xest;

            Y.row(i) = y.transpose();
        }
        return Y;
    }
};

int main(int argc, char** argv) {
    Kalman kalman;
    MatrixXd x = MatrixXd::Random(4, 2);

    auto estimated = kalman.estimate(x);

    cout << "x=\n" << x << endl;
    cout << "estimated=\n" << estimated << endl;
    return 0;
}
