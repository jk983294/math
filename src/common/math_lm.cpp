#include <math_lm.h>
#include <math_random.h>
#include <math_stats.h>
#include <math_vector.h>
#include <cmath>
#include <cstdio>
#include <stats.hpp>
#include <tuple>
#include <Eigen/Dense>
#include <zerg_template.h>

namespace ornate {
namespace {
constexpr int SKIP_INDEX_NAN = -1;
constexpr int SKIP_INDEX_TEST = 0;
constexpr int SKIP_INDEX_TRAIN = 1;

static std::vector<double> calc_p_value(const Eigen::MatrixXd& A, const Eigen::VectorXd& residual, const Eigen::VectorXd& coef) {
    const int n = A.rows();
    const int k = A.cols();
    double s2 = residual.dot(residual) / (n - k);
    auto& pinv = (A.transpose() * A).completeOrthogonalDecomposition().pseudoInverse();
    const Eigen::VectorXd& std_err = (s2 * pinv.diagonal()).array().sqrt();
    const Eigen::VectorXd& tvalue = coef.array() / std_err.array();
    std::vector<double> pvalue(tvalue.data(), tvalue.data() + tvalue.size());
    for (size_t j = 0; j < pvalue.size(); ++j) {
        if (std::isfinite(pvalue[j])) {
            pvalue[j] = 2 * (1 - stats::pt(std::abs(pvalue[j]), n - tvalue.size()));
        } else {
            printf("pvalue inf found %zu, %f\n", j, pvalue[j]);
        }
    }
    return pvalue;
}

static double calc_r2(const Eigen::MatrixXd& X_mat, const Eigen::VectorXd& b,
                           const Eigen::VectorXd& coef) {
    Eigen::VectorXd ytot = b.array() - b.array().mean();
    Eigen::VectorXd fitted = X_mat * coef;
    Eigen::VectorXd residual = b - fitted;
    return 1 - residual.dot(residual) / ytot.dot(ytot);
}
}  // namespace
double LmModel::train_lm_first_n_features(int feature_count) {
    printf("train step 1\n");
    selected.clear();
    for (int i = 0; i < feature_count; ++i) selected.push_back(i);

    size_t train_cnt, test_cnt;
    std::tie(train_cnt, test_cnt) = calc_skip_indices();

    m_coefs.clear();
    m_signal_mean = {0};  // 截距项
    m_signal_sd = {0};
    int col = feature_count;
    if (m_param.m_intercept) col += 1;
    Eigen::MatrixXd X(train_cnt, col);
    Eigen::MatrixXd X_test(test_cnt, col);
    double* pX = X.data();
    double* pX_test = X_test.data();
    if (m_param.m_intercept) {
        std::fill(pX, pX + train_cnt, 1.0);           // intercept
        std::fill(pX_test, pX_test + test_cnt, 1.0);  // intercept
    }

    printf("train step 2\n");
    Eigen::VectorXd b(train_cnt);
    Eigen::VectorXd b_test(test_cnt);
    fill_data(train_cnt, SKIP_INDEX_TRAIN, m_param.m_y, b.data());
    fill_data(test_cnt, SKIP_INDEX_TEST, m_param.m_y, b_test.data());
    Eigen::VectorXd ytot = b.array() - b.array().mean();

    int X_offset = 0;
    if (m_param.m_intercept) X_offset = 1;
    for (int i = 0; i < feature_count; ++i) {
        double mean_ = ornate::mean(m_param.m_features[i], m_param.n_row);
        double sd_ = ornate::std(m_param.m_features[i], m_param.n_row);
        m_signal_mean.push_back(mean_);
        m_signal_sd.push_back(sd_);

        fill_data(train_cnt, SKIP_INDEX_TRAIN, m_param.m_features[i], &pX[(i + X_offset) * train_cnt]);
        fill_data(test_cnt, SKIP_INDEX_TEST, m_param.m_features[i], &pX[(i + X_offset) * test_cnt]);
    }

    // LOG_INFO("train_whole start");
    Eigen::VectorXd coef = X.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
    const Eigen::VectorXd& fitted = X * coef;
    const Eigen::VectorXd& residual = b - fitted;
    double r_squared = 1 - residual.dot(residual) / ytot.dot(ytot);

    std::vector<double> pvalue = calc_p_value(X, residual, coef);

    printf("train step 3\n");
    if (m_param.m_intercept) {
        m_coefs.resize(coef.size());
        std::copy(coef.data(), coef.data() + coef.size(), m_coefs.begin());
    } else {
        m_coefs.resize(col + 1);
        m_coefs[0] = 0.;  // intercept
        std::copy(coef.data(), coef.data() + col, m_coefs.begin() + 1);
    }

    double r_squared_test = NAN;
    if (m_test_cnt > 0) {
        r_squared_test = calc_r2(X_test, b_test, coef);
    }

    printf("train_whole r_squared=(%f,%f)\n", r_squared, r_squared_test);
    printf("f: intercept,%s\n", ztool::head(m_param.m_f_names, 0).c_str());
    printf("coef: %s\n", ztool::head(m_coefs, 0).c_str());
    printf("pvalue: %s\n", ztool::head(pvalue, 0).c_str());
    return r_squared;
}

std::pair<size_t, size_t> LmModel::calc_skip_indices() {
    m_skip_indices.clear();
    m_skip_indices.resize(m_param.n_row, SKIP_INDEX_NAN);
    auto untradable = *m_param.m_untradable;
    ornate::MyRandom myRandom;
    size_t feature_all_nan_count{0}, feature_any_nan_count{0}, na_cnt{0}, train_cnt{0}, test_cnt{0}, untradable_cnt{0};
    for (size_t i = 0; i < m_param.n_row; ++i) {
        if (!untradable.empty() && untradable.size() == m_param.n_row && untradable[i]) {
            ++untradable_cnt;
        } else if (!std::isfinite(m_param.m_y[i])) {
            ++na_cnt;
        } else if (!if_any_feature_valid(i)) {
            ++feature_all_nan_count;
        } else if (m_param.m_skip_na_feature && !if_all_feature_valid(i)) {
            ++feature_any_nan_count;
        } else {
            if (need_test_set(i)) {
                m_skip_indices[i] = SKIP_INDEX_TEST;
                ++test_cnt;
            } else {
                m_skip_indices[i] = SKIP_INDEX_TRAIN;
                ++train_cnt;
            }
        }
    }
    printf("na=%zu, untradable=%zu, feature_nan=(%zu,%zu), train/test=(%zu,%zu), total=%zu\n",
           na_cnt, untradable_cnt, feature_all_nan_count, feature_any_nan_count, train_cnt, test_cnt, m_param.n_row);
    return {train_cnt, test_cnt};
}

void LmModel::train(const TrainParam& param) {
    m_param = param;
    if (m_param.is_step) {
        train_lm_step();
    } else {
        train_lm_whole();
    }
}

std::vector<double> LmModel::fit_new(size_t n, const std::unordered_map<std::string, const double*>& name2features) {
    // TODO
    return {};
}

double LmModel::fit_row(const std::vector<double>& inputs) {
    int f_cnt = m_param.m_features.size();
    double ret = m_coefs[0];
    for (int i = 0; i < f_cnt; ++i) {
        if (ornate::isvalid(inputs[i]))
            ret += inputs[i] * m_coefs[i + 1];
    }
    return ret;
}

std::vector<double> LmModel::fit() {
    std::vector<double> rets(m_param.n_row, NAN);
    int f_cnt = m_param.m_features.size();
    std::vector<double> inputs(f_cnt, NAN);
    for (size_t i = 0; i < m_param.n_row; ++i) {
        for (int j = 0; j < f_cnt; ++j) {
            inputs[j] = m_param.m_features[j][i];
        }
        rets[i] = fit_row(inputs);
    }
    return rets;
}

void LmModel::train_lm_step() {}

double LmModel::train_lm_whole() {
    if (m_max_feature_num > 0 && m_max_feature_num < (int)m_param.m_features.size()) {
        return train_lm_first_n_features(m_max_feature_num);
    } else {
        return train_lm_first_n_features(m_param.m_features.size());
    }
}

/**
 * 根据y来填x, 如果y invalid，该行不要
 * 如果x invalid, 取0
 */
void LmModel::fill_data(size_t expected_row_cnt, int flag, const double* src, double* dest) {
    size_t offset = 0;
    for (size_t i = 0; i < m_param.n_row; ++i) {
        if (m_skip_indices[i] == flag) {
            if (ornate::isvalid(src[i]))
                dest[offset] = src[i];
            else
                dest[offset] = 0;

            ++offset;
        }
    }
    if (offset != expected_row_cnt) {
        printf("na cnt wrong %zu %zu\n", offset, expected_row_cnt);
        exit(1);
    }
}

bool LmModel::if_feature_valid(size_t row_id) {
    int f_count = selected.size();
    for (int j = 0; j < f_count; ++j) {
        int fid = selected[j];
        if (std::isfinite((m_param.m_features[fid])[row_id])) return true;
    }
    return false;
}

bool LmModel::if_any_feature_valid(size_t row_id) {
    int f_count = m_param.m_features.size();
    for (int fid = 0; fid < f_count; ++fid) {
        if (std::isfinite((m_param.m_features[fid])[row_id])) return true;
    }
    return false;
}
bool LmModel::if_all_feature_valid(size_t row_id) {
    int f_count = m_param.m_features.size();
    for (int fid = 0; fid < f_count; ++fid) {
        if (!std::isfinite((m_param.m_features[fid])[row_id])) return false;
    }
    return true;
}
}  // namespace ornate