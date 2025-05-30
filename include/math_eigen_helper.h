#pragma once

#include <Eigen/Dense>
#include <vector>

namespace ornate {
inline Eigen::MatrixXd cov2corr(const Eigen::MatrixXd& cov);
inline Eigen::MatrixXd corr2cov(const Eigen::MatrixXd& corr, const Eigen::VectorXd& sd);
inline std::vector<double> ToVector(const Eigen::VectorXd& vec);
inline std::vector<double> ToVector(const Eigen::MatrixXd& m);
inline Eigen::VectorXd ToVector(const std::vector<double>& vec);
inline Eigen::MatrixXd ToMat(const std::vector<double>& vec, size_t nrow, size_t ncol, bool row_major, bool symmetric);
inline void append(Eigen::VectorXd& target, double val);
inline void append(Eigen::VectorXd& target, const Eigen::VectorXd& vec);
inline void append(Eigen::MatrixXd& target, const Eigen::MatrixXd& m, bool vertically = true);
inline void append(Eigen::MatrixXd& target, const Eigen::VectorXd& vec, bool vertically = true);
inline void resize(Eigen::VectorXd& target, size_t new_size);
inline void resize(Eigen::MatrixXd& target, size_t new_row, size_t new_col);

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

inline Eigen::MatrixXd cov2corr(const Eigen::MatrixXd& cov_) {
    int n = cov_.rows();

    Eigen::MatrixXd corr_ = Eigen::MatrixXd::Zero(n, n);
    Eigen::VectorXd sd = cov_.diagonal().array().sqrt();

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            corr_(i, j) = cov_(i, j) / (sd(i) * sd(j));
        }
    }
    return corr_;
}

inline Eigen::MatrixXd corr2cov(const Eigen::MatrixXd& corrMatrix, const Eigen::VectorXd& sd) {
    int n = corrMatrix.rows();

    Eigen::MatrixXd cov_ = Eigen::MatrixXd::Zero(n, n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cov_(i, j) = corrMatrix(i, j) * sd(i) * sd(j);
        }
    }
    return cov_;
}

inline Eigen::VectorXd ToVector(const std::vector<double>& vec) {
    Eigen::VectorXd _vector = Eigen::Map<const Eigen::VectorXd>(vec.data(), vec.size());
    return _vector;
}

inline std::vector<double> ToVector(const Eigen::VectorXd& vec) {
    std::vector<double> ret(vec.size());
    const double* ptr = vec.data();
    std::copy(ptr, ptr + vec.size(), ret.data());
    return ret;
}

inline std::vector<double> ToVector(const Eigen::MatrixXd& m) {
    std::vector<double> ret(m.rows() * m.cols());
    const double* ptr = m.data();
    std::copy(ptr, ptr + ret.size(), ret.data());
    return ret;
}

inline Eigen::MatrixXd ToMat(const std::vector<double>& vec, size_t nrow, size_t ncol, bool row_major, bool symmetric) {
    symmetric = (nrow == ncol) & symmetric;
    if (row_major && !symmetric) {
        Eigen::Map<const Eigen::MatrixXd> m(vec.data(), ncol, nrow);
        return m.transpose();
    } else {
        return Eigen::Map<const Eigen::MatrixXd>(vec.data(), nrow, ncol);
    }
}

inline void append(Eigen::VectorXd& target, double val) {
    target.conservativeResize(target.size() + 1);
    target(target.size() - 1) = val;
}

inline void append(Eigen::VectorXd& target, const Eigen::VectorXd& vec) {
    if (target.size() == 0) {
        target = vec;
    } else {
        target.conservativeResize(target.size() + vec.size());
        target.tail(vec.size()) = vec;
    }
}

inline void append(Eigen::MatrixXd& target, const Eigen::MatrixXd& m, bool vertically) {
    if (vertically) {
        if (target.rows() == 0) {
            target = m;
        } else {
            target.conservativeResize(target.rows() + m.rows(), target.cols());
            target.bottomRows(m.rows()) = m;
        }
    } else {
        if (target.cols() == 0) {
            target = m;
        } else {
            target.conservativeResize(target.rows(), target.cols() + m.cols());
            target.rightCols(m.cols()) = m;
        }
    }
}

inline void append(Eigen::MatrixXd& target, const Eigen::VectorXd& vec, bool vertically) {
    if (vertically) {
        if (target.rows() == 0) {
            target.conservativeResize(1, vec.size());
        } else {
            target.conservativeResize(target.rows() + 1, target.cols());
        }
        target.row(target.rows() - 1) = vec.transpose();
    } else {
        if (target.rows() == 0) {
            target.conservativeResize(vec.size(), 1);
        } else {
            target.conservativeResize(target.rows(), target.cols() + 1);
        }
        target.col(target.cols() - 1) = vec;
    }
}

inline void resize(Eigen::VectorXd& target, size_t new_size) {
    Eigen::VectorXd resized_(new_size);
    resized_.setZero();
    resized_.head(target.size()) = target;
    target.swap(resized_);
}

inline void resize(Eigen::MatrixXd& target, size_t new_row, size_t new_col) {
    Eigen::MatrixXd resized_(new_row, new_col);
    resized_.setZero();
    resized_.block(0, 0, target.rows(), target.cols()) = target;
    target.swap(resized_);
}
}  // namespace ornate
