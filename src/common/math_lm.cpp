#include <math_lm.h>
#include <math_random.h>
#include <math_stats.h>
#include <math_vector.h>
#include <armadillo>
#include <cmath>
#include <cstdio>
#include <stats.hpp>
#include <tuple>

namespace ornate {
namespace {
constexpr int SKIP_INDEX_NAN = -1;
constexpr int SKIP_INDEX_TEST = 0;
constexpr int SKIP_INDEX_TRAIN = 1;

static std::vector<double> calc_p_value(const arma::mat& A, const arma::colvec& residual, const arma::vec& coef) {
    double s2 = arma::dot(residual, residual) / int(A.n_rows - A.n_cols);
    const arma::colvec& std_err = arma::sqrt(s2 * arma::diagvec(arma::pinv(arma::trans(A) * A)));
    const arma::colvec& tvalue = coef / std_err;
    std::vector<double> pvalue(tvalue.begin(), tvalue.end());
    for (size_t j = 0; j < pvalue.size(); ++j) {
        if (std::isfinite(pvalue[j])) {
            pvalue[j] = 2 * (1 - stats::pt(std::abs(pvalue[j]), A.n_rows - tvalue.size()));
        } else {
            printf("pvalue inf found %zu, %f\n", j, pvalue[j]);
        }
    }
    return pvalue;
}

static double calc_r2(std::vector<double>& X, int row, int col, const arma::colvec& b, const arma::colvec& ytot,
                      arma::vec& coef1) {
    const arma::mat A1(X.data(), row, col, false, false);
    bool status1 = solve(coef1, A1, b);
    const arma::colvec& fitted1 = A1 * coef1;
    const arma::colvec& residual1 = b - fitted1;
    return status1 ? 1 - arma::dot(residual1, residual1) / arma::dot(ytot, ytot) : NAN;
}

static double calc_r2_test(std::vector<double>& X_test, int row, int col, const arma::colvec& b_test,
                           const arma::vec& coef) {
    const arma::colvec& ytot_test = b_test - arma::mean(b_test);
    const arma::mat A_test(X_test.data(), row, col, false, false);
    const arma::colvec& fitted_test = A_test * coef;
    const arma::colvec& residual_test = b_test - fitted_test;
    return 1 - arma::dot(residual_test, residual_test) / arma::dot(ytot_test, ytot_test);
}
}  // namespace
double LmModel::train_lm_first_n_features(int feature_count) {
    selected.clear();
    for (int i = 0; i < feature_count; ++i) selected.push_back(i);

    size_t train_cnt, test_cnt;
    std::tie(train_cnt, test_cnt) = calc_skip_indices();

    m_coefs.clear();
    m_signal_mean = {0};  // 截距项
    m_signal_sd = {0};
    int col = feature_count;
    if (m_intercept) col += 1;
    vector<double> X(col * train_cnt);
    vector<double> X_test(col * test_cnt);
    if (m_intercept) {
        std::fill(X.begin(), X.begin() + train_cnt, 1.0);           // intercept
        std::fill(X_test.begin(), X_test.begin() + test_cnt, 1.0);  // intercept
    }

    vector<double> real_y(train_cnt), real_y_test(test_cnt);
    fill_data(train_cnt, SKIP_INDEX_TRAIN, m_y, real_y.data());
    fill_data(test_cnt, SKIP_INDEX_TEST, m_y, real_y_test.data());
    const arma::colvec b(real_y.data(), real_y.size(), false, false);
    const arma::colvec b_test(real_y_test.data(), real_y_test.size(), false, false);
    const arma::colvec& ytot = b - arma::mean(b);

    int X_offset = 0;
    if (m_intercept) X_offset = 1;
    for (int i = 0; i < feature_count; ++i) {
        double mean_ = ornate::mean(m_features[i], total_row);
        double sd_ = ornate::std(m_features[i], total_row);
        m_signal_mean.push_back(mean_);
        m_signal_sd.push_back(sd_);

        fill_data(train_cnt, SKIP_INDEX_TRAIN, m_features[i], &X[(i + X_offset) * train_cnt]);
        fill_data(test_cnt, SKIP_INDEX_TEST, m_features[i], &X[(i + X_offset) * test_cnt]);
    }

    // LOG_INFO("train_whole start");
    const arma::mat A(X.data(), train_cnt, col, false, false);
    arma::vec coef;
    bool status = solve(coef, A, b);
    if (not status) {
        printf("solve failed status=%d\n", status);
        return NAN;
    }
    const arma::colvec& fitted = A * coef;
    const arma::colvec& residual = b - fitted;
    double r_squared = 1 - arma::dot(residual, residual) / arma::dot(ytot, ytot);

    if (m_compress_with_p_value) {
        std::vector<double> pvalue = calc_p_value(A, residual, coef);
        if (compress_with_p_value(pvalue, X, train_cnt, X_test, test_cnt, m_intercept)) {
            arma::vec coef1;
            col = (int)selected.size();
            if (m_intercept) col += 1;
            double r_squared1 = calc_r2(X, train_cnt, col, b, ytot, coef1);
            if (std::isfinite(r_squared1)) {
                printf("after t value compress r2=%f, before=%f\n", r_squared1, r_squared);
            } else {
                printf("solve failed after t value compress\n");
            }
            coef = coef1;
        }
    }

    if (m_intercept) {
        m_coefs.resize(col);
        std::copy(coef.begin(), coef.end(), m_coefs.begin());
    } else {
        m_coefs.resize(col + 1);
        m_coefs[0] = 0.;  // intercept
        std::copy(coef.begin(), coef.end(), m_coefs.begin() + 1);
    }

    double r_squared_test = NAN;
    if (m_test_cnt > 0) {
        r_squared_test = calc_r2_test(X_test, test_cnt, col, b_test, coef);
    }

    printf("train_whole status=%d, r_squared=(%f,%f)\n", status, r_squared, r_squared_test);
    return r_squared;
}

bool LmModel::compress_with_p_value(const std::vector<double>& pvalue, vector<double>& X, int n, vector<double>& X_test,
                                    int n_test, bool has_intercept) {
    vector<int> selected_tmp;
    std::vector<double> signal_mean_tmp = {0}, signal_sd_tmp = {0};

    for (size_t i = has_intercept ? 1 : 0; i < pvalue.size(); ++i) {
        int fid = selected[i - (has_intercept ? 1 : 0)];
        if (pvalue[i] < m_p_value_threshold) {
            selected_tmp.push_back(fid);
            signal_mean_tmp.push_back(m_signal_mean[i + (has_intercept ? 0 : 1)]);
            signal_sd_tmp.push_back(m_signal_sd[i + (has_intercept ? 0 : 1)]);
            auto to_offset = (selected_tmp.size() - 1 + (has_intercept ? 1 : 0));
            auto from_offset = i;
            if (to_offset != from_offset) {
                double* to_ = X.data() + n * to_offset;
                double* from_ = X.data() + n * from_offset;
                std::copy(from_, from_ + n, to_);

                double* to_test = X_test.data() + n_test * to_offset;
                double* from_test = X_test.data() + n_test * from_offset;
                std::copy(from_test, from_test + n_test, to_test);
            }
        } else {
            printf("discard feature %d, t-value=%f\n", fid, pvalue[i]);
        }
    }

    printf("compress, before=%zu, after=%zu\n", selected.size(), selected_tmp.size());
    bool compressed = false;
    if (selected_tmp.size() < selected.size()) {
        compressed = true;
    }
    selected.swap(selected_tmp);
    m_signal_mean.swap(signal_mean_tmp);
    m_signal_sd.swap(signal_sd_tmp);
    return compressed;
}

std::pair<size_t, size_t> LmModel::calc_skip_indices() {
    m_skip_indices.clear();
    m_skip_indices.resize(total_row, SKIP_INDEX_NAN);
    ornate::MyRandom myRandom;
    size_t feature_all_nan_count = 0, feature_any_nan_count = 0, na_cnt = 0, train_cnt = 0, test_cnt = 0;
    for (size_t i = 0; i < total_row; ++i) {
        if (!std::isfinite(m_y[i])) {
            ++na_cnt;
        } else if (!if_any_feature_valid(i)) {
            ++feature_all_nan_count;
        } else if (m_skip_na_feature && !if_all_feature_valid(i)) {
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
    printf("na=%zu, feature_nan=(%zu,%zu), train/test=(%zu,%zu), total=%zu\n", na_cnt, feature_all_nan_count,
           feature_any_nan_count, train_cnt, test_cnt, total_row);
    return {train_cnt, test_cnt};
}

void LmModel::lm(bool is_step) {
    if (is_step) {
        train_lm_step();
    } else {
        train_lm_whole();
    }
}

void LmModel::train_lm_step() {}

double LmModel::train_lm_whole() {
    if (m_max_feature_num > 0 && m_max_feature_num < (int)m_features.size()) {
        return train_lm_first_n_features(m_max_feature_num);
    } else {
        return train_lm_first_n_features(m_features.size());
    }
}

/**
 * 根据y来填x, 如果y invalid，该行不要
 * 如果x invalid, 取0
 */
void LmModel::fill_data(size_t expected_row_cnt, int flag, const double* src, double* dest) {
    size_t offset = 0;
    for (size_t i = 0; i < total_row; ++i) {
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
        if (std::isfinite((m_features[fid])[row_id])) return true;
    }
    return false;
}

bool LmModel::if_any_feature_valid(size_t row_id) {
    int f_count = m_features.size();
    for (int fid = 0; fid < f_count; ++fid) {
        if (std::isfinite((m_features[fid])[row_id])) return true;
    }
    return false;
}
bool LmModel::if_all_feature_valid(size_t row_id) {
    int f_count = m_features.size();
    for (int fid = 0; fid < f_count; ++fid) {
        if (!std::isfinite((m_features[fid])[row_id])) return false;
    }
    return true;
}
}  // namespace ornate