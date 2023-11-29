#include <math_lm.h>
#include <math_random.h>
#include <math_vector.h>
#include <cmath>
#include <cstdio>
#include <stats.hpp>
#include <tuple>
#include <Eigen/Dense>
#include <zerg_template.h>
#include <fstream>
#include <zerg_file.h>
#include <math_lasso.h>

namespace ornate {
namespace {
constexpr int SKIP_INDEX_NAN = -1;
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
}  // namespace
double LmModel::train_lm() {
    Eigen::MatrixXd X;
    Eigen::VectorXd y;
    int col = prepare_Xy(X, y);

    // LOG_INFO("train_whole start");
    Eigen::VectorXd coef = X.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(y);
    const Eigen::VectorXd& fitted = X * coef;
    const Eigen::VectorXd& residual = y - fitted;
    Eigen::VectorXd ytot = y.array() - y.array().mean();
    double r_squared = 1 - residual.dot(residual) / ytot.dot(ytot);

    m_pvalues = calc_p_value(X, residual, coef);

    if (m_param.m_intercept) {
        m_coefs.resize(coef.size());
        std::copy(coef.data(), coef.data() + coef.size(), m_coefs.begin());
    } else {
        m_coefs.resize(col + 1);
        m_coefs[0] = 0.;  // intercept
        std::copy(coef.data(), coef.data() + col, m_coefs.begin() + 1);
    }

    printf("train_whole r_squared=%f\n", r_squared);
    printf("f: intercept,%s\n", ztool::head(m_param.m_f_names, 0).c_str());
    printf("coef: %s\n", ztool::head(m_coefs, 0).c_str());
    printf("pvalue: %s\n", ztool::head(m_pvalues, 0).c_str());
    return r_squared;
}

int LmModel::prepare_Xy(Eigen::MatrixXd& X, Eigen::VectorXd& y) {
    int feature_count = m_param.m_features.size();
    size_t train_cnt = calc_skip_indices();

    m_coefs.clear();
    int col = feature_count;
    if (m_param.m_intercept) col += 1;
    X.resize(train_cnt, col);
    double* pX = X.data();
    if (m_param.m_intercept) {
        std::fill(pX, pX + train_cnt, 1.0);           // intercept
    }

    y.resize(train_cnt);
    fill_data(train_cnt, SKIP_INDEX_TRAIN, m_param.m_y, y.data());


    int X_offset = 0;
    if (m_param.m_intercept) X_offset = 1;
    for (int i = 0; i < feature_count; ++i) {
        fill_data(train_cnt, SKIP_INDEX_TRAIN, m_param.m_features[i], &pX[(i + X_offset) * train_cnt]);
    }
    return col;
}

double LmModel::train_lasso() {
    Eigen::MatrixXd X;
    Eigen::VectorXd y;
    int col = prepare_Xy(X, y);

    // LOG_INFO("train_whole start");
    ornate::LassoModel model(m_param.lasso_lambda, m_param.lasso_n_iter, m_param.lasso_error_threshold,
                             m_param.m_lasso_need_scale, false);
    model.set_verbose(m_param.m_verbose);
    model.fit(X, y);
    Eigen::VectorXd y_hat = model.predict(X);
    Eigen::VectorXd coef = model.get_model().coef;
    Eigen::VectorXd ytot = y.array() - y.array().mean();

    const Eigen::VectorXd& fitted = X * coef;
    const Eigen::VectorXd& residual = y - fitted;
    double r_squared = 1 - residual.dot(residual) / ytot.dot(ytot);

    m_pvalues = calc_p_value(X, residual, coef);

    if (m_param.m_intercept) {
        m_coefs.resize(coef.size());
        std::copy(coef.data(), coef.data() + coef.size(), m_coefs.begin());
    } else {
        m_coefs.resize(col + 1);
        m_coefs[0] = 0.;  // intercept
        std::copy(coef.data(), coef.data() + col, m_coefs.begin() + 1);
    }

    printf("train_whole r_squared=%f\n", r_squared);
    printf("f: intercept,%s\n", ztool::head(m_param.m_f_names, 0).c_str());
    printf("coef: %s\n", ztool::head(m_coefs, 0).c_str());
    printf("pvalue: %s\n", ztool::head(m_pvalues, 0).c_str());
    return r_squared;
}

size_t LmModel::calc_skip_indices() {
    m_skip_indices.clear();
    m_skip_indices.resize(m_param.n_row, SKIP_INDEX_NAN);
    auto untradable = *m_param.m_untradable;
    ornate::MyRandom myRandom;
    size_t feature_all_nan_count{0}, feature_any_nan_count{0}, na_cnt{0}, train_cnt{0}, untradable_cnt{0};
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
            m_skip_indices[i] = SKIP_INDEX_TRAIN;
            ++train_cnt;
        }
    }
    printf("na=%zu, untradable=%zu, feature_nan=(%zu,%zu), train=%zu, total=%zu\n",
           na_cnt, untradable_cnt, feature_all_nan_count, feature_any_nan_count, train_cnt, m_param.n_row);
    return train_cnt;
}

void LmModel::train(const TrainParam& param) {
    m_param = param;
    if (param.m_lasso) {
        train_lasso();
    } else {
        train_lm();
    }
}

std::vector<double> LmModel::fit_new(size_t n, const std::unordered_map<std::string, const double*>& name2features) {
    std::vector<double> rets(n, NAN);
    int f_cnt = m_param.m_f_names.size();
    std::vector<const double*> _features(f_cnt, nullptr);
    for (int j = 0; j < f_cnt; ++j) {
        auto itr = name2features.find(m_param.m_f_names[j]);
        if (itr == name2features.end()) {
            printf("fit_new failed missing %s\n", m_param.m_f_names[j].c_str());
            return rets;
        } else {
            _features[j] = itr->second;
        }
    }
    std::vector<double> inputs(f_cnt, NAN);
    for (size_t i = 0; i < n; ++i) {
        for (int j = 0; j < f_cnt; ++j) {
            inputs[j] = _features[j][i];
        }
        rets[i] = fit_row(inputs);
    }
    return rets;
}

double LmModel::fit_row(const std::vector<double>& inputs) {
    int f_cnt = inputs.size();
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

bool LmModel::save(std::string path) {
    path = ztool::FileExpandUser(path);
    std::ofstream ofs(path, std::ofstream::out | std::ofstream::trunc);

    if (!ofs) {
        printf("save lm model open file %s failed\n", path.c_str());
        return false;
    }
    if (m_param.m_f_names.empty()) return false;
    if (m_param.m_lasso) ofs << "#LASSO\n";
    else ofs << "#LM\n";
    ofs << "F_CNT," << m_param.m_f_names.size() << "\n";
    ofs << "INTERCEPT," << m_coefs.front() << "\n";
    for (size_t i = 0; i < m_param.m_f_names.size(); ++i) {
        ofs << m_param.m_f_names[i] << "," << m_coefs[i + 1] << "\n";
    }
    printf("write lm model to %s\n", path.c_str());
    return true;
}
bool LmModel::load(std::string path) {
    path = ztool::FileExpandUser(path);
    std::ifstream ifs(path, std::ifstream::in);

    if (!ifs) {
        printf("load lm model open file %s failed\n", path.c_str());
        return false;
    }

    string s;
    int f_cnt = 0;
    while (getline(ifs, s)) {
        if (s.empty() || s.front() == '#') continue;
        auto itr = s.find_first_of(',');
        auto first = s.substr(0, itr);
        auto second = s.substr(itr + 1);
        if (first == "F_CNT") {
            f_cnt = std::stoi(second);
        } else if (first == "INTERCEPT") {
            m_coefs.push_back(std::stod(second));
        } else {
            m_coefs.push_back(std::stod(second));
            m_param.m_f_names.push_back(first);
        }
    }

    ifs.close();
    return (int)m_param.m_f_names.size() == f_cnt;
}
}  // namespace ornate