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

std::vector<double> Model::fit_new(TrainType type, size_t n, const std::unordered_map<std::string, const double*>& name2features) {
    if (type == TrainType::LM)
        return m_lm.fit_new(n, name2features);
    else {
        return {};
    }
}

std::vector<double> Model::fit(TrainType type) {
    if (type == TrainType::LM)
        return m_lm.fit();
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

    if (type == TrainType::LM)
        m_lm.train(m_param);
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

}  // namespace ornate