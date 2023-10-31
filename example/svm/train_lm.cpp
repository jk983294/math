#include <cstring>
#include <random>
#include <math_lm.h>

using namespace ornate;


int main(int argc, char** argv) {
    Model model;
    size_t n = 100;
    size_t total_row = n * n;
    std::random_device rd;
    std::mt19937 generator(rd());
    std::normal_distribution<double> nd(0., 0.001);
    model.total_row = total_row;

    std::vector<double> x0(total_row, NAN);
    std::vector<double> x1(total_row, NAN);
    std::vector<double> y(total_row, NAN);

    // y = 2 * x0 + 3 * x1  - 1 + delta
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            x0[i * n + j] = i;
            x1[i * n + j] = j;
            y[i * n + j] = 2. * i + 3. * j - 1 + nd(generator);
        }
    }

    model.m_y = y.data();
    model.m_features.push_back(x0.data());
    model.m_features.push_back(x1.data());
    model.fit_lm();
    return 0;
}
