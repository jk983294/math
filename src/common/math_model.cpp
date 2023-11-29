#include <math_lm.h>
#include <math_vector.h>
#include <cstdio>
#include <tuple>

namespace ornate {

std::vector<double> Model::fit_new(TrainType type, size_t n, const std::unordered_map<std::string, const double*>& name2features) {
    if (type == TrainType::LM)
        return m_lm.fit_new(n, name2features);
    else if (type == TrainType::SVM)
        return fit_new_svm(n, name2features);
    else {
        return {};
    }
}

std::vector<double> Model::fit(TrainType type) {
    if (type == TrainType::LM || type == TrainType::Lasso)
        return m_lm.fit();
    else if (type == TrainType::SVM)
        return fit_svm();
    else {
        return {};
    }
}

bool Model::train(TrainType type, const std::vector<std::string>& features) {
    m_param.m_features.clear();
    m_param.m_f_names.clear();
    m_param.m_untradable = &m_untradable;

    if (features.empty()) {
        m_param.m_f_names = m_f_names;
        m_param.m_features = m_features;
    } else {
        for (auto& name : features) {
            auto itr = m_name2features.find(name);
            if (itr == m_name2features.end() || itr->second == nullptr) {
                printf("WARN: cannot find feature %s\n", name.c_str());
                return false;
            } else {
                m_param.m_features.push_back(itr->second);
                m_param.m_f_names.push_back(itr->first);
            }
        }
    }

    if (m_param.m_features.empty() || m_param.n_row == 0) {
        printf("WARN: no feature/row to train train\n");
        return false;
    }
    m_param.m_y = m_y;

    if (type == TrainType::LM) {
        m_param.m_lasso = false;
        m_lm.train(m_param);
    } else if (type == TrainType::SVM)
        train_svm();
    else if (type == TrainType::Lasso) {
        m_param.m_lasso = true;
        m_lm.train(m_param);
    }
    else {
        printf("WARN: un support train type %d\n", type);
        return false;
    }
    return true;
}

Model::~Model() {
    for (auto* f : m_real_datum) delete f;
    m_real_datum.clear();
}

std::vector<double> Model::fit_new_svm(size_t n, const std::unordered_map<std::string, const double*>& name2features) {
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
        rets[i] = m_svm.predict_x(inputs.data());
    }
    return rets;
}

std::vector<double> Model::fit_svm() {
    auto untradable = *m_param.m_untradable;
    std::vector<double> predicts(m_param.n_row, NAN);
    std::vector<double> xs(m_param.m_features.size());
    for (size_t i = 0; i < m_param.n_row; ++i) {
        double y = m_param.m_y[i];
        if (not std::isfinite(y)) continue;
        if (not untradable.empty() && untradable[i]) continue;

        for (size_t j = 0; j < m_param.m_features.size(); ++j) {
            double x = m_param.m_features[j][i];
            if (not std::isfinite(x)) x = 0;
            xs[j] = x;
        }

        predicts[i] = m_svm.predict_x(xs.data());
    }
    return predicts;
}

void Model::train_svm() {
    auto untradable = *m_param.m_untradable;
    size_t valid_n = 0;
    for (size_t i = 0; i < m_param.n_row; ++i) {
        double y = m_param.m_y[i];
        if (not std::isfinite(y)) continue;
        if (not untradable.empty() && untradable[i]) continue;
        valid_n++;
    }

    ornate::DataSet ds;
    ds.reserve(valid_n, m_param.m_features.size());
    std::vector<double> xs(m_param.m_features.size());
    for (size_t i = 0; i < m_param.n_row; ++i) {
        double y = m_param.m_y[i];
        if (not std::isfinite(y)) continue;
        if (not untradable.empty() && untradable[i]) continue;
        for (size_t j = 0; j < m_param.m_features.size(); ++j) {
            double x = m_param.m_features[j][i];
            if (not std::isfinite(x)) x = 0;
            xs[j] = x;
        }
        ds.add_row(y, xs);
    }

    m_svm.param.solver_type = 11;
    m_svm.flag_solver_specified = true;
    m_svm.flag_find_parameters = true;
    m_svm.bias = -1; // no bias term
    m_svm.init();
    m_svm.set_data(ds);
    m_svm.work();
}

bool Model::add_feature_real(std::string name, std::vector<double> f) {
    auto* f_real = new std::vector<double>;
    f_real->swap(f);
    m_real_datum.push_back(f_real);
    return add_feature(name, f_real->data());
}
bool Model::add_feature_real(std::string name, std::vector<bool> f) {
    auto* f_real = new std::vector<bool>;
    f_real->swap(f);
    m_real_bool_datum.push_back(f_real);
    return add_feature(name, f_real);
}
void Model::set_y_real(std::vector<double> y) {
    auto* y_real = new std::vector<double>;
    y_real->swap(y);
    m_real_datum.push_back(y_real);
    set_y(y_real->data());
}
bool Model::add_feature(std::string name, const std::vector<bool>* f) {
    auto itr = m_name2bFeatures.find(name);
    if (itr == m_name2bFeatures.end()) {
        m_name2bFeatures[name] = f;
        m_f_bNames.push_back(name);
        m_bFeatures.push_back(f);
        return true;
    } else {
        printf("WARN: dupe bool feature %s\n", name.c_str());
        return false;
    }
}
bool Model::add_feature(std::string name, const double* f) {
    auto itr = m_name2features.find(name);
    if (itr == m_name2features.end()) {
        m_name2features[name] = f;
        m_f_names.push_back(name);
        m_features.push_back(f);
        return true;
    } else {
        printf("WARN: dupe double feature %s\n", name.c_str());
        return false;
    }
}

void Model::reset_untradable() {
    m_untradable.resize(m_param.n_row);
    std::fill(m_untradable.begin(), m_untradable.end(), false);
}
void Model::set_untradable(const std::vector<bool>& untradable) {
    if (m_untradable.empty()) reset_untradable();

    for (size_t i = 0; i < m_param.n_row; ++i) {
        if (untradable[i]) m_untradable[i] = true;
    }
}
void Model::set_untradable(std::string name) {
    auto itr = m_name2bFeatures.find(name);
    if (itr == m_name2bFeatures.end()) {
        set_untradable(*itr->second);
    } else {
        printf("WARN: untradable feature %s not found\n", name.c_str());
    }
}

bool Model::save(TrainType type, std::string path) {
    if (type == TrainType::LM || type == TrainType::Lasso)
        return m_lm.save(path);
    else if (type == TrainType::SVM)
        return m_svm.save(path, m_param.m_f_names);
    else {
        printf("Model::save type not support\n");
        return false;
    }
}

bool Model::load(TrainType type, std::string path) {
    if (type == TrainType::LM || type == TrainType::Lasso)
        return m_lm.load(path);
    else if (type == TrainType::SVM)
        return m_svm.load(path, m_param.m_f_names);
    else {
        printf("Model::load type not support\n");
        return false;
    }
}
}  // namespace ornate