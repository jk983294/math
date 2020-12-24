#include <math_random.h>
#include <algorithm>
#include <armadillo>
#include <cassert>
#include <cstdio>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;
using namespace arma;

enum AssessCriterion { AdjustR2, AIC };

class ForwardStep {
public:
    void add_data(int id, double* data) {
        if (m_id2data.find(id) == m_id2data.end()) {
            m_id2data[id] = data;
        } else {
            printf("add_data failed, id %d exist", id);
        }
    }
    void set_y(double* data) { y = data; }

    void process() {
        purge_data();
        m_selected.clear();
        m_candidates.clear();
        for (const auto& item : m_id2data) m_candidates.insert(item.first);
        printf("cnt=%zu\n", m_candidates.size());

        for (const auto& id : m_candidates) {
            double v = calc_criterion(id);
            printf("id=%d score=%f\n", id, v);
        }
    }

private:
    void purge_data() {
        purge_data(y);
        for (auto& item : m_id2data) purge_data(item.second);
    }
    void purge_data(double* data) {
        for (size_t i = 0; i < n; ++i) {
            if (!std::isfinite(data[i])) data[i] = 0.;
        }
    }

    const double* get_data_pointer(int id) {
        auto* ptr = m_id2data[id];
        assert(ptr != nullptr);
        return ptr;
    }
    double calc_criterion(int id) {
        std::size_t p = m_selected.size() + 1;
        size_t col_cnt = p;
        if (has_intercept) ++col_cnt;

        mat A(n, col_cnt);
        double* pA = A.memptr();

        size_t col_idx = 0;
        for (; col_idx < m_selected.size(); ++col_idx) {
            auto* pData = get_data_pointer(m_selected[col_idx]);
            std::copy(pData, pData + n, pA + col_idx * n);
        }

        {
            auto* pData = get_data_pointer(id);
            std::copy(pData, pData + n, pA + col_idx * n);
            ++col_idx;
        }

        if (has_intercept) {
            std::fill_n(pA + col_idx * n, n, 1.0);
        }
        vec b(y, n);
        vec coef;
        A.print("A:");
        bool status = solve(coef, A, b);
        if (!status) return NAN;

        const arma::colvec& fitted = A * coef;
        const arma::colvec& residual = b - fitted;

        if (m_criterion == AdjustR2) {
            double r_squared;
            if (has_intercept) {
                const arma::colvec& ytot = b - arma::mean(b);
                r_squared = 1 - arma::dot(residual, residual) / arma::dot(ytot, ytot);
            } else {
                r_squared = 1 - arma::dot(residual, residual) / arma::dot(b, b);
            }

            double adjusted_r_squared = 1 - (1 - r_squared) * double(n - 1) / double(n - p - 1);
            return adjusted_r_squared;
        } else {  // AIC = 2k + n * ln(RSS)
            double aic = 2.0 * p + n * std::log(arma::dot(residual, residual));
            return -1.0 * aic;  // 越大越好, 和 adjusted_r_squared 方向一致
        }
        // return NAN;
    }

public:
    double* y{nullptr};
    std::vector<int> m_selected;
    std::unordered_set<int> m_candidates;
    std::unordered_map<int, double*> m_id2data;
    size_t n{0};  // observation count
    bool has_intercept{false};
    AssessCriterion m_criterion{AssessCriterion::AdjustR2};
};

int main() {
    size_t n = 5;
    ForwardStep fs;
    // fs.has_intercept = true;
    fs.n = n;
    vector<double> ret = ornate::generate_uniform_float<double>(n, 0.0, 1.0);
    vector<double> f0 = ornate::generate_uniform_float<double>(n, 0.0, 1.0);
    vector<double> f1 = ornate::generate_uniform_float<double>(n, 0.0, 1.0);
    //    vector<double> f2 = ornate::generate_uniform_float<double>(n, 0.0, 1.0);
    //    vector<double> f3 = ornate::generate_uniform_float<double>(n, 0.0, 1.0);
    //    vector<double> f4 = ornate::generate_uniform_float<double>(n, 0.0, 1.0);
    //    vector<double> f5 = ornate::generate_uniform_float<double>(n, 0.0, 1.0);
    fs.add_data(0, f0.data());
    fs.add_data(0, f1.data());
    //    fs.add_data(2, f2.data());
    //    fs.add_data(3, f3.data());
    //    fs.add_data(4, f4.data());
    //    fs.add_data(5, f5.data());
    fs.set_y(ret.data());
    fs.process();
    return 0;
}
