#pragma once

#include <Eigen/Dense>

namespace ornate {
inline Eigen::VectorXd eigen_col_wise_std(const Eigen::MatrixXd& a) {
    Eigen::VectorXd mean = a.colwise().sum() / a.rows();
    return ((a.rowwise() - mean.transpose()).array().pow(2).colwise().sum() / (a.rows() - 1)).sqrt();
}
inline void eigen_col_norm(Eigen::MatrixXd& a, const Eigen::VectorXd& stddev) {
    a.array().rowwise() /= stddev.transpose().array();
}
inline Eigen::MatrixXd eigen_col_norm_copy(const Eigen::MatrixXd& a, const Eigen::VectorXd& stddev) {
    return a.array().rowwise() / stddev.transpose().array();
}

inline Eigen::VectorXd eigen_row_wise_std(const Eigen::MatrixXd& a) {
    Eigen::VectorXd mean = a.rowwise().sum() / a.cols();
    return ((a.colwise() - mean).array().pow(2).rowwise().sum() / (a.cols() - 1)).sqrt();
}
inline void eigen_row_norm(Eigen::MatrixXd& a, const Eigen::VectorXd& stddev) { a.array().colwise() /= stddev.array(); }
inline Eigen::MatrixXd eigen_row_norm_copy(const Eigen::MatrixXd& a, const Eigen::VectorXd& stddev) {
    return a.array().colwise() / stddev.array();
}
inline void eigen_vec_replace(Eigen::VectorXd& x, double from, double to, double eps = 1e-6) {
    double* pd = x.data();
    for (long i = 0; i < x.size(); ++i) {
        if (std::abs(pd[i] - from) < eps) {
            pd[i] = to;
        }
    }
}
/**
 * @param a
 * @param dim 0 row 1 column
 */
inline Eigen::VectorXd eigen_mean(const Eigen::MatrixXd& a, int dim) {
    if (dim == 0)
        return a.rowwise().mean();
    else
        return a.colwise().mean();
}

inline double eigen_dot(const Eigen::VectorXd& a, const Eigen::VectorXd& b) { return a.transpose() * b; }

template <typename T>
inline auto eigen_ele_mult(const T& a, const T& b) {
    return a.array() * b.array();
}

inline void eigen_remove_row(Eigen::MatrixXd& matrix, unsigned int rowToRemove) {
    unsigned int numRows = matrix.rows() - 1;
    unsigned int numCols = matrix.cols();

    if (rowToRemove < numRows)
        matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) = matrix.bottomRows(numRows - rowToRemove);

    matrix.conservativeResize(numRows, numCols);
}

inline void eigen_remove_column(Eigen::MatrixXd& matrix, unsigned int colToRemove) {
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols() - 1;

    if (colToRemove < numCols)
        matrix.block(0, colToRemove, numRows, numCols - colToRemove) = matrix.rightCols(numCols - colToRemove);

    matrix.conservativeResize(numRows, numCols);
}

inline double calc_r2(const Eigen::MatrixXd& X_mat, const Eigen::VectorXd& y, const Eigen::VectorXd& coef) {
    Eigen::VectorXd ytot = y.array() - y.array().mean();
    Eigen::VectorXd fitted = X_mat * coef;
    Eigen::VectorXd residual = y - fitted;
    return 1 - residual.dot(residual) / ytot.dot(ytot);
}
}  // namespace ornate
