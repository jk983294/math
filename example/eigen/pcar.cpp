#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

using namespace Eigen;
using namespace std;


ArrayXXd cumsum(const ArrayXXd& arr) {
    ArrayXXd result(arr.rows(), arr.cols());
    result.setZero();
    for (int i = 0; i < arr.rows(); ++i) {
        double sum_ = 0;
        for (int j = 0; j < arr.cols(); ++j) {
            sum_ += arr(i, j);
            result(i, j) = sum_;
        }
    }
    return result;
}

pair<int, MatrixXd> pca(const MatrixXd& X, double variance_threshold) {
    // Compute the covariance matrix
    MatrixXd cov = X.adjoint() * X / (X.rows() - 1);

    // Compute the eigenvalues and eigenvectors
    SelfAdjointEigenSolver<MatrixXd> eigensolver(cov);
    MatrixXd eigenvalues = eigensolver.eigenvalues();
    MatrixXd eigenvectors = eigensolver.eigenvectors();

    cout << "eigenvalues: " << eigenvalues << endl;

    // Sort eigenvalues and eigenvectors in descending order
    ArrayXXd sorted_eigenvalues = eigenvalues.transpose().array().reverse();
    MatrixXd sorted_eigenvectors = eigenvectors.rowwise().reverse();

    cout << sorted_eigenvalues << std::endl;

    ArrayXXd cumulative_variance = cumsum(sorted_eigenvalues) / sorted_eigenvalues.sum();

    cout << "cumulative_variance: " << cumulative_variance << std::endl;

    int num_components = 0;
    while (num_components < cumulative_variance.size() && cumulative_variance(num_components) < variance_threshold) {
        num_components++;
    }
    num_components++; // Include the component that crosses the threshold

    // Select the top 'num_components' principal components
    MatrixXd pca_loading = sorted_eigenvectors.leftCols(num_components);

    cout << "num_components: " << num_components << std::endl;
    cout << "pca_loading:\n" << pca_loading << std::endl;

    return {num_components, pca_loading};
}

// Function to perform linear regression
VectorXd linear_regression(const MatrixXd& X, const VectorXd& y) {
    MatrixXd X_pinv = (X.adjoint() * X).inverse() * X.adjoint();
    VectorXd beta = X_pinv * y;
    return beta;
}

int main() {
    int n_obv = 100;
    int n_variable = 10;
    MatrixXd X = MatrixXd::Random(n_obv, n_variable);
    VectorXd y = VectorXd::Random(n_obv);

    auto [num_components, pca_loading] = pca(X, 0.5);

    // Transform the data to the new space
    MatrixXd transformed_data = X * pca_loading;

    // Perform linear regression using the transformed data
    VectorXd beta = linear_regression(transformed_data, y);

    // Interpret the results in terms of the original variables
    MatrixXd beta_original = pca_loading * beta;

    cout << "Regression coefficients in terms of original variables:\n" << beta_original << endl;

    return 0;
}