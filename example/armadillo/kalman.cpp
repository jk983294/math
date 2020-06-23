#include <armadillo>
#include <iostream>

using namespace std;
using namespace arma;

class Kalman {
private:
    mat A;  // state transition matrix
    mat H;  // measurement matrix
    mat Q, R;
    mat xest, pest;  // state conditions
    double dt;

public:
    // constructor, sets up data structures
    Kalman() : dt(1.0) {
        A.eye(6, 6);
        A(0, 2) = A(1, 3) = A(2, 4) = A(3, 5) = dt;
        H.zeros(2, 6);
        H(0, 0) = H(1, 1) = 1.0;

        Q.eye(6, 6);
        R = 1000 * eye(2, 2);

        xest.zeros(6, 1);
        pest.zeros(6, 6);
    }

    // sole member function: estimate model
    mat estimate(const mat& Z) {
        size_t n = Z.n_rows, k = Z.n_cols;
        mat Y = zeros(n, k);
        mat xprd, pprd, S, B, kalmangain;
        colvec y;

        for (unsigned int i = 0; i < n; i++) {
            colvec z = Z.row(i).t();
            // predicted state and covariance
            xprd = A * xest;
            pprd = A * pest * A.t() + Q;

            // estimation
            S = H * pprd.t() * H.t() + R;
            B = H * pprd.t();

            kalmangain = trans(solve(S, B));

            // estimated state and covariance
            xest = xprd + kalmangain * (z - H * xprd);
            pest = pprd - kalmangain * H * pprd;

            // compute the estimated measurements
            y = H * xest;

            Y.row(i) = y.t();
        }
        return Y;
    }
};

int main(int argc, char** argv) {
    Kalman kalman;
    mat x = randu<mat>(4, 2);

    auto estimated = kalman.estimate(x);

    cout << "x=\n" << x << endl;
    cout << "estimated=\n" << estimated << endl;
    return 0;
}
