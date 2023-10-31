#pragma once

#include <vector>

using std::vector;

namespace ornate {
struct LmModel {
    void lm(bool is_step);
    std::vector<double*> m_features;
    double* m_y{nullptr};
    std::vector<double*> m_real;
    size_t total_row;
    bool m_intercept{true};
    bool m_skip_na_feature{false};  // default only skip y and all x invalid, enable this will skip any x invalid
    bool m_compress_with_p_value{false};
    double m_p_value_threshold{0.1};
    vector<int> selected;
    std::vector<int> m_skip_indices;
    int m_max_feature_num{-1};
    size_t m_train_cnt{0}, m_test_offset{0}, m_test_cnt{0};
    std::vector<double> m_coefs;
    std::vector<double> m_signal_mean, m_signal_sd;

private:
    void train_lm_step();
    double train_lm_whole();
    double train_lm_first_n_features(int feature_count);

private:
    bool if_feature_valid(size_t row_id);
    bool if_any_feature_valid(size_t row_id);
    bool if_all_feature_valid(size_t row_id);
    std::pair<size_t, size_t> calc_skip_indices();
    bool need_test_set(size_t idx) const { return m_test_cnt > 0 && idx >= m_test_offset; }
    bool compress_with_p_value(const std::vector<double>& pvalue, vector<double>& X, int n, vector<double>& X_test,
                               int n_test, bool has_intercept);
    void fill_data(size_t expected_row_cnt, int flag, const double* src, double* dest);
};

struct Model {
    void fit_lm(bool intercept = true);

    std::vector<double*> m_features;
    double* m_y{nullptr};
    size_t total_row;
    LmModel m_lm;
};
}  // namespace ornate
