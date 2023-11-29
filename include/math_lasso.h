#pragma once

#include <math_eigen_helper.h>

using namespace std;

namespace ornate {
using Eigen::MatrixXd;
using Eigen::RowVectorXd;
using Eigen::VectorXd;

struct LModel {
    RowVectorXd mean;
    RowVectorXd sd;
    VectorXd coef;
};

class LassoModel {
private:
    class StandardScaler {
    public:
        void fit(const MatrixXd &X) {
            m_mean_ = X.colwise().mean();
            auto t = X.rowwise() - m_mean_;
            m_sd_ = (t.array() * t.array()).matrix().colwise().sum().array().sqrt();
        }
        MatrixXd transform(const MatrixXd &X) {  // x - mean(x) / sd(x)
            MatrixXd t(X);
            for (int i = 0; i < X.cols(); ++i) t.col(i) /= m_sd_(i);
            return t;
        }
        MatrixXd fit_transform(const MatrixXd &X) {
            fit(X);
            return transform(X);
        }
        void set_param(const RowVectorXd &mean, const RowVectorXd &var) {
            m_mean_ = mean;
            m_sd_ = var;
        }
        RowVectorXd get_mean() { return m_mean_; }
        RowVectorXd get_sd() { return m_sd_; }

    private:
        RowVectorXd m_mean_;
        RowVectorXd m_sd_;
    };

public:
    LassoModel(double lambda, size_t n_iter, double e, bool need_scale = true, bool need_bias = true)
        : lambda_(lambda), n_iter_(n_iter), e_(e), m_need_scale{need_scale}, m_need_bias{need_bias} {
        if (m_need_scale) scalaer = new StandardScaler();
    }
    ~LassoModel() { delete scalaer; }
    VectorXd predict(const MatrixXd &X) { return transform_x(X) * coef_; }
    void fit(const MatrixXd &X, const VectorXd &y) {
        MatrixXd X_train = transform_x(X);
        coef_ = VectorXd::Zero(X_train.cols());
        for (size_t iter = 0; iter < n_iter_; ++iter) {
            RowVectorXd z = (X_train.array() * X_train.array()).colwise().sum();
            VectorXd tmp = VectorXd::Zero(X_train.cols());
            assert(z.rows() == tmp.cols());
            for (auto k = 0; k < X_train.cols(); ++k) {
                double wk = coef_(k);
                coef_(k) = 0;
                double p_k = X_train.col(k).transpose() * (y - X_train * coef_);
                double w_k = 0.0;
                if (p_k < -lambda_ / 2.0)
                    w_k = (p_k + lambda_ / 2.0) / z(k);
                else if (p_k > lambda_ / 2.0)
                    w_k = (p_k - lambda_ / 2.0) / z(k);
                else
                    w_k = 0.0;
                tmp(k) = w_k;
                coef_(k) = wk;
            }
            double erroe_ = (coef_ - tmp).norm();
            if (m_verbose) {
                printf("round=%zu, e=%f\n", iter, erroe_);
            }
            if (erroe_ < e_) {
                if (m_verbose) {
                    printf("early stop round %zu, threshold=%f, e=%f\n", iter, e_, erroe_);
                }
                break;
            }

            coef_ = tmp;
        }
    }

    void set_model(const VectorXd &coef, const RowVectorXd &mean, const RowVectorXd &var) {
        coef_ = coef;
        if (scalaer) scalaer->set_param(mean, var);
    }
    LModel get_model() {
        LModel model;
        model.coef = coef_;
        if (scalaer && m_need_scale) {
            model.mean = scalaer->get_mean();
            model.sd = scalaer->get_sd();
        }
        return model;
    }

    void set_verbose(bool flag) { m_verbose = flag; }

private:
    MatrixXd transform_x(const MatrixXd &X) {
        MatrixXd X_test = X;
        if (m_need_scale) {
            X_test = scalaer->fit_transform(X);
        }
        if (m_need_bias) {
            X_test = combine_bias(X_test);
        }
        return X_test;
    }
    MatrixXd combine_bias(const MatrixXd &X) {
        MatrixXd t(X.rows(), X.cols() + 1);
        t.col(X.cols()) = VectorXd::Ones(X.rows());
        t.leftCols(X.cols()) = X;
        return t;
    }

private:
    double lambda_;
    size_t n_iter_;
    double e_;
    StandardScaler *scalaer{nullptr};
    VectorXd coef_;
    bool m_verbose{false};
    bool m_need_scale{true};
    bool m_need_bias{true};
};
}  // namespace ornate
