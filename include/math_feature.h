#ifndef ORNATE_MATH_FEATURE_H
#define ORNATE_MATH_FEATURE_H

#include <armadillo>
#include <cassert>
#include <cstdio>
#include <string>
#include <unordered_map>
#include <unordered_set>

namespace ornate {
class ForwardStep {
public:
    enum AssessCriterion { AdjustR2, AIC };

    ~ForwardStep() {
        for(auto* container_ : m_containers) delete container_;
        m_containers.clear();
    }
    void add_data(int id, double* data) {
        if (n == 0) {
            printf("must set y first along with n\n");
            return;
        }

        if (m_id2data.find(id) == m_id2data.end()) {
            if (remove_na_ret) {
                auto& container_ = get_container();
                size_t na_idx = 0;
                for (size_t i = 0; i < original_n; ++i) {
                    if (na_idx < m_na_index.size() && i == m_na_index[na_idx]) {
                        ++na_idx;
                    } else {
                        // purge
                        if(std::isfinite(data[i])) {
                            container_.push_back(data[i]);
                        } else {
                            container_.push_back(0.0);
                        }
                    }
                }
                m_id2data[id] = container_.data();
            } else {
                m_id2data[id] = data;
            }
        } else {
            printf("add_data failed, id %d exist\n", id);
        }
    }
    void set_y(double* data, size_t n_) {
        original_n = n_;
        if (remove_na_ret) {
            for (size_t i = 0; i < n_; ++i) {
                if(std::isfinite(data[i])) {
                    y_container.push_back(data[i]);
                } else {
                    m_na_index.push_back(i);
                }
            }
            y = y_container.data();
            n = y_container.size();
        } else {
            y = data;
            n = n_;
        }
    }

    void process() {
        if (not remove_na_ret) {
            purge_data();
        }
        m_selected.clear();
        m_candidates.clear();
        for (const auto& item : m_id2data) m_candidates.push_back(item.first);
        printf("cnt=%zu\n", m_candidates.size());

        while (!m_candidates.empty()) {
            double curr_best = NAN;
            int curr_best_id = std::numeric_limits<int>::max();
            int choice_cnt = (int)m_candidates.size();
            std::vector<double> scores(choice_cnt);

            #pragma omp parallel for schedule(dynamic, 1) num_threads(thread_num)
            for (int i = 0; i < choice_cnt; ++i) {
                int id = m_candidates[i];
                scores[i] = calc_criterion(id);
            }

            for (int i = 0; i < choice_cnt; ++i) {
                // printf("id=%d, score=%f\n", m_candidates[i], scores[i]);
                if (std::isnan(curr_best) || curr_best < scores[i]) {
                    curr_best = scores[i];
                    curr_best_id = m_candidates[i];
                }
            }

            printf("current round id=%d, best=%f\n", curr_best_id, curr_best);
            if (std::isnan(history_best) || history_best < curr_best) {
                auto itr = std::find(m_candidates.begin(), m_candidates.end(), curr_best_id);
                if(itr != m_candidates.end()) {
                    m_selected.push_back(curr_best_id);
                    m_candidates.erase(itr);
                    history_best = curr_best;
                } else {
                    break;
                }
            } else {
                break;
            }
        }
    }

    void set_thread_num(int num) { thread_num = num; }
    void set_remove_na_ret(bool remove_na_ret_) { remove_na_ret = remove_na_ret_; }
    void set_has_intercept(bool has_intercept_) { has_intercept = has_intercept_; }
    void set_criterion(const std::string& criterion) {
        if (criterion == "aic") m_criterion = AssessCriterion::AIC;
        else m_criterion = AssessCriterion::AdjustR2;
    }
    std::vector<int> get_selected() { return m_selected; }

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

    /**
     * @return score 越大越好, 和 adjusted_r_squared 方向一致
     */
    double calc_criterion(int id) {
        std::size_t p = m_selected.size() + 1;
        size_t col_cnt = p;
        if (has_intercept) ++col_cnt;

        arma::mat A(n, col_cnt);
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
        arma::vec b(y, n);
        arma::vec coef;
        // A.print("A:");
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
    }

    std::vector<double>& get_container() {
        m_containers.push_back(new std::vector<double>());
        return *m_containers.back();
    }

private:
    double* y{nullptr};
    std::vector<int> m_selected;
    std::vector<int> m_candidates;
    std::unordered_map<int, double*> m_id2data;
    size_t original_n{0}; // from ret vector
    size_t n{0};  // observation count
    int thread_num{1};
    bool has_intercept{false};
    bool remove_na_ret{false};
    double history_best{NAN};
    std::vector<std::size_t> m_na_index;
    std::vector<double> y_container;
    std::vector<std::vector<double>*> m_containers;
    AssessCriterion m_criterion{AssessCriterion::AdjustR2};
};
}

#endif
