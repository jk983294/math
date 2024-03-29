#include <math_lm.h>
#include <math_stats.h>
#include <cstring>
#include <random>

using namespace ornate;

void train_lm(Model& model, std::vector<double>& y, std::vector<double>& x0) {
    model.train(TrainType::LM, {});
    auto fitted = model.fit(TrainType::LM);
    double pcor0 = ornate::corr(fitted, y);
    model.train(TrainType::LM, {"x0"});
    fitted = model.fit(TrainType::LM);
    double pcor1 = ornate::corr(fitted, y);
    printf("corr %f <-> %f\n", pcor0, pcor1);

    std::string m_path = "/tmp/test.lm.model";
    model.save(TrainType::LM, m_path);

    Model model1;
    model1.load(TrainType::LM, m_path);
    std::unordered_map<std::string, const double*> name2features;
    name2features["x0"] = x0.data();
    auto fitted1 = model1.fit_new(TrainType::LM, x0.size(), name2features);
    double pcor2 = ornate::corr(fitted1, y);
    printf("corr %f <-> %f\n", pcor0, pcor2);
}

void train_svm(Model& model, std::vector<double>& y, std::vector<double>& x0) {
    model.train(TrainType::SVM, {});
    auto fitted = model.fit(TrainType::SVM);
    double pcor0 = ornate::corr(fitted, y);
    model.train(TrainType::SVM, {"x0"});
    fitted = model.fit(TrainType::SVM);
    double pcor1 = ornate::corr(fitted, y);
    printf("corr %f <-> %f\n", pcor0, pcor1);

    std::string m_path = "/tmp/test.svm.model";
    model.save(TrainType::SVM, m_path);

    Model model1;
    model1.load(TrainType::SVM, m_path);
    std::unordered_map<std::string, const double*> name2features;
    name2features["x0"] = x0.data();
    auto fitted1 = model1.fit_new(TrainType::SVM, x0.size(), name2features);
    double pcor2 = ornate::corr(fitted1, y);
    printf("corr %f <-> %f\n", pcor0, pcor2);
}

void train_lasso(Model& model, std::vector<double>& y, std::vector<double>& x0) {
    model.m_param.m_lasso_need_scale = false;
    model.train(TrainType::Lasso, {});
    auto fitted = model.fit(TrainType::Lasso);
    double pcor0 = ornate::corr(fitted, y);
    model.train(TrainType::Lasso, {"x0"});
    fitted = model.fit(TrainType::Lasso);
    double pcor1 = ornate::corr(fitted, y);
    printf("corr %f <-> %f\n", pcor0, pcor1);

    std::string m_path = "/tmp/test.lasso.model";
    model.save(TrainType::Lasso, m_path);

    Model model1;
    model1.load(TrainType::Lasso, m_path);
    std::unordered_map<std::string, const double*> name2features;
    name2features["x0"] = x0.data();
    auto fitted1 = model1.fit_new(TrainType::Lasso, x0.size(), name2features);
    double pcor2 = ornate::corr(fitted1, y);
    printf("corr %f <-> %f\n", pcor0, pcor2);
}

int main(int argc, char** argv) {
    Model model;
    size_t n = 100;
    size_t n_row = n * n;
    std::random_device rd;
    std::mt19937 generator(rd());
    std::normal_distribution<double> nd(0., 0.001);

    std::vector<double> x0(n_row, NAN);
    std::vector<double> x1(n_row, NAN);
    std::vector<double> y(n_row, NAN);

    // y = 2 * x0 + 3 * x1  - 1 + delta
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            x0[i * n + j] = i;
            x1[i * n + j] = j;
            y[i * n + j] = 2. * i + 3. * j - 1 + nd(generator);
        }
    }

    model.set_n(n_row);
    //    model.set_y(y.data());
    //    model.add_feature("x0", x0.data());
    //    model.add_feature("x1", x1.data());
    model.set_y_real(y);
    model.add_feature_real("x0", x0);
    model.add_feature_real("x1", x1);
    vector<bool> untradable(n_row, false);
    model.set_untradable(untradable);
    // train_lm(model, y, x0);
    // train_svm(model, y, x0);
    train_lasso(model, y, x0);
    return 0;
}
