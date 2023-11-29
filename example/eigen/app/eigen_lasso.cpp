#include <math_lasso.h>
#include <math_stats.h>
#include <iostream>

using namespace std;

int main(int argc, char* argv[]) {
    std::random_device rd;
    std::mt19937 generator(rd());
    std::normal_distribution<double> nd(0., 0.001);

    // y = 2 * x0 + 3 * x1 + delta
    size_t n = 1000;
    Eigen::MatrixXd matX(n * n, 3);
    Eigen::VectorXd y(n * n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            double x0 = i, x1 = j;
            double y_val = 2 * x0 + 3 * x1 + nd(generator);
            y(i * n + j) = y_val;
            matX(i * n + j, 0) = x0;
            matX(i * n + j, 1) = x1;
            matX(i * n + j, 2) = nd(generator) * 100000;
        }
    }

    if (y.size() != matX.rows()) cout << "Number of responses must be equal to number of rows of X!" << endl;

    ornate::LassoModel model(0.2, 100, 0.01, false, false);
    model.set_verbose(true);
    model.fit(matX, y);
    Eigen::VectorXd y_hat = model.predict(matX);
    Eigen::VectorXd beta = model.get_model().coef;
    cout << "beta=" << beta << endl;

    double pcor = ornate::corr(y_hat.data(), y.data(), n * n);
    printf("corr %f\n", pcor);
}
