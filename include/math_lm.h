#pragma once

#include <linear_svm.h>
#include <math_lasso.h>
#include <cmath>
#include <unordered_map>
#include <vector>

using std::vector;

namespace ornate {

enum TrainType : int { LM = 0, Lasso, SVM };

struct TrainParam {
    std::vector<const double*> m_features;
    std::vector<std::string> m_f_names;
    const double* m_y{nullptr};
    const std::vector<bool>* m_untradable{nullptr};
    size_t n_row{0};

    bool m_intercept{true};
    bool m_skip_na_feature{false};  // if any feature nan will cause this row skipped
    bool m_verbose{true};
    bool m_lasso{false};
    bool m_lasso_need_scale{false};
    // bool m_lasso_need_bias{false}; same with m_intercept
    double lasso_lambda{0.2};
    double lasso_error_threshold{0.01};
    size_t lasso_n_iter{100};

    svm::parameter m_svm_param;
};

struct LmModel {
    void train(const TrainParam& param);
    std::vector<double> fit_new(size_t n, const std::unordered_map<std::string, const double*>& name2features);
    std::vector<double> fit();
    bool save(std::string path);
    bool load(std::string path);

    TrainParam m_param;
    std::vector<int> m_skip_indices;
    size_t m_train_cnt{0};
    std::vector<double> m_coefs, m_pvalues;

private:
    double fit_row(const std::vector<double>& inputs);
    double train_lm();
    double train_lasso();

private:
    bool if_any_feature_valid(size_t row_id);
    bool if_all_feature_valid(size_t row_id);
    size_t calc_skip_indices();
    void fill_data(size_t expected_row_cnt, int flag, const double* src, double* dest);
    int prepare_Xy(Eigen::MatrixXd& X, Eigen::VectorXd& y);
};

struct Model {
    ~Model();
    bool train(TrainType type, const std::vector<std::string>& features);
    std::vector<double> fit_new(TrainType type, size_t n,
                                const std::unordered_map<std::string, const double*>& name2features);
    std::vector<double> fit(TrainType type);

    bool add_feature(std::string name, const double* f);
    bool add_feature(std::string name, const std::vector<bool>* f);
    void set_y(const double* y) { m_y = y; }
    void set_n(size_t n) { m_param.n_row = n; }
    size_t get_n() const { return m_param.n_row; }

    bool add_feature_real(std::string name, std::vector<double> f);
    bool add_feature_real(std::string name, std::vector<bool> f);
    void set_y_real(std::vector<double> y);

    void set_intercept(bool flag) { m_param.m_intercept = flag; }
    bool get_intercept() const { return m_param.m_intercept; }
    void set_skip_na_feature(bool flag) { m_param.m_skip_na_feature = flag; }
    bool get_skip_na_feature() const { return m_param.m_skip_na_feature; }
    std::vector<std::string> get_feature_names() const { return m_f_names; }

    void reset_untradable();
    void set_untradable(const std::vector<bool>& untradable);
    void set_untradable(std::string name);
    LmModel& get_lm_model() { return m_lm; }

    bool save(TrainType type, std::string path);
    bool load(TrainType type, std::string path);

private:
    void train_svm();
    std::vector<double> fit_svm();
    std::vector<double> fit_new_svm(size_t n, const std::unordered_map<std::string, const double*>& name2features);

public:
    std::vector<std::vector<double>*> m_real_datum;
    std::vector<std::vector<bool>*> m_real_bool_datum;
    std::vector<bool> m_untradable;
    std::vector<const double*> m_features;
    std::vector<std::string> m_f_names;
    std::vector<const std::vector<bool>*> m_bFeatures;
    std::vector<std::string> m_f_bNames;
    std::unordered_map<std::string, const double*> m_name2features;
    std::unordered_map<std::string, const std::vector<bool>*> m_name2bFeatures;
    const double* m_y{nullptr};
    LmModel m_lm;
    svm::LinearSvm m_svm;
    TrainParam m_param;
};
}  // namespace ornate
