#include "linear_svm.h"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <fstream>

#ifdef __cplusplus
extern "C" {
#endif

extern double dnrm2_(int *, double *, int *);
extern double ddot_(int *, double *, int *, double *, int *);
extern int daxpy_(int *, double *, double *, int *, double *, int *);
extern int dscal_(int *, double *, double *, int *);

#ifdef __cplusplus
}
#endif

namespace ornate::svm {

class sparse_operator {
public:
    static double nrm2_sq(const double *x, int n) {
        double ret = 0;
        for (int i = 0; i < n; ++i) {
            ret += x[i] * x[i];
        }
        return ret;
    }
    static double dot(const double *s, const double *x, int n) {
        double ret = 0;
        for (int i = 0; i < n; ++i) {
            ret += s[i] * x[i];
        }
        return ret;
    }
    static void axpy(const double a, const double *x, double *y, int n) {
        for (int i = 0; i < n; ++i) {
            y[i] += a * x[i];
        }
    }
};

// L2-regularized empirical risk minimization
// min_w w^Tw/2 + \sum C_i \xi(w^Tx_i), where \xi() is the loss

class l2r_erm_fun : public function {
public:
    l2r_erm_fun(const problem *prob, const parameter *param, double *C);
    ~l2r_erm_fun() override = default;

    double fun(double *w) override;
    double linesearch_and_update(double *w, double *s, double *f, double *g, double alpha) override;
    int get_nr_variable() override;

protected:
    virtual double C_times_loss(int i, double wx_i) = 0;
    void Xv(double *v, double *Xv);

    double *C{nullptr};
    const problem *prob{nullptr};
    std::vector<double> wx;   // size = #y
    std::vector<double> tmp;  // size = #y, a working array
    double wTw{0};
    int regularize_bias;
};

l2r_erm_fun::l2r_erm_fun(const problem *prob, const parameter *param, double *C) {
    this->prob = prob;

    wx.resize(prob->l);
    tmp.resize(prob->l);
    this->C = C;
    this->regularize_bias = param->regularize_bias;
}

double l2r_erm_fun::fun(double *w) {
    int i;
    double f = 0;
    int l = prob->l;
    int w_size = get_nr_variable();

    wTw = 0;
    Xv(w, wx.data());

    for (i = 0; i < w_size; i++) wTw += w[i] * w[i];
    if (regularize_bias == 0) wTw -= w[w_size - 1] * w[w_size - 1];
    for (i = 0; i < l; i++) f += C_times_loss(i, wx[i]);
    f = f + 0.5 * wTw;

    return f;
}

int l2r_erm_fun::get_nr_variable() { return prob->n; }

// On entry *f must be the function value of w
// On exit w is updated and *f is the new function value
double l2r_erm_fun::linesearch_and_update(double *w, double *s, double *f, double *g, double alpha) {
    int i;
    int l = prob->l;
    double sTs = 0;
    double wTs = 0;
    double gTs = 0;
    double eta = 0.01;
    int w_size = get_nr_variable();
    int max_num_linesearch = 20;
    double fold = *f;
    Xv(s, tmp.data());

    for (i = 0; i < w_size; i++) {
        sTs += s[i] * s[i];
        wTs += s[i] * w[i];
        gTs += s[i] * g[i];
    }
    if (regularize_bias == 0) {
        // bias not used in calculating (w + \alpha s)^T (w + \alpha s)
        sTs -= s[w_size - 1] * s[w_size - 1];
        wTs -= s[w_size - 1] * w[w_size - 1];
    }

    int num_linesearch = 0;
    for (num_linesearch = 0; num_linesearch < max_num_linesearch; num_linesearch++) {
        double loss = 0;
        for (i = 0; i < l; i++) {
            double inner_product = tmp[i] * alpha + wx[i];
            loss += C_times_loss(i, inner_product);
        }
        *f = loss + (alpha * alpha * sTs + wTw) / 2.0 + alpha * wTs;
        if (*f - fold <= eta * alpha * gTs) {
            for (i = 0; i < l; i++) wx[i] += alpha * tmp[i];
            break;
        } else
            alpha *= 0.5;
    }

    if (num_linesearch >= max_num_linesearch) {
        *f = fold;
        return 0;
    } else
        for (i = 0; i < w_size; i++) w[i] += alpha * s[i];

    wTw += alpha * alpha * sTs + 2 * alpha * wTs;
    return alpha;
}

void l2r_erm_fun::Xv(double *v, double *Xv) {
    for (int i = 0; i < prob->l; i++) {
        Xv[i] = sparse_operator::dot(v, prob->m_x[i], prob->n);
    }
}

class l2r_l2_svc_fun : public l2r_erm_fun {
public:
    l2r_l2_svc_fun(const problem *prob, const parameter *param, double *C);
    ~l2r_l2_svc_fun() override = default;

    void grad(double *w, double *g) override;
    void Hv(double *s, double *Hs) override;

    void get_diag_preconditioner(double *M) override;

protected:
    void subXTv(double *v, double *XTv);
    std::vector<int> I;  // #y
    int sizeI{0};

private:
    double C_times_loss(int i, double wx_i) override;
};

l2r_l2_svc_fun::l2r_l2_svc_fun(const problem *prob, const parameter *param, double *C) : l2r_erm_fun(prob, param, C) {
    I.resize(prob->l);
}

double l2r_l2_svc_fun::C_times_loss(int i, double wx_i) {
    double d = 1 - prob->m_y[i] * wx_i;
    if (d > 0)
        return C[i] * d * d;
    else
        return 0;
}

void l2r_l2_svc_fun::grad(double *w, double *g) {
    int i;
    auto &y = prob->m_y;
    int l = prob->l;
    int w_size = get_nr_variable();

    sizeI = 0;
    for (i = 0; i < l; i++) {
        tmp[i] = wx[i] * y[i];
        if (tmp[i] < 1) {
            tmp[sizeI] = C[i] * y[i] * (tmp[i] - 1);
            I[sizeI] = i;
            sizeI++;
        }
    }
    subXTv(tmp.data(), g);

    for (i = 0; i < w_size; i++) g[i] = w[i] + 2 * g[i];
    if (regularize_bias == 0) g[w_size - 1] -= w[w_size - 1];
}

void l2r_l2_svc_fun::get_diag_preconditioner(double *M) {
    int i;
    int w_size = get_nr_variable();

    for (i = 0; i < w_size; i++) M[i] = 1;
    if (regularize_bias == 0) M[w_size - 1] = 0;

    for (i = 0; i < sizeI; i++) {
        int idx = I[i];
        auto *xi = prob->m_x[idx];
        for (int j = 0; j < w_size; ++j) {
            M[j] += xi[j] * xi[j] * C[idx] * 2;
        }
    }
}

void l2r_l2_svc_fun::Hv(double *s, double *Hs) {
    int i;
    int w_size = get_nr_variable();

    for (i = 0; i < w_size; i++) Hs[i] = 0;
    for (i = 0; i < sizeI; i++) {
        auto xi = prob->m_x[I[i]];
        double xTs = sparse_operator::dot(s, xi, w_size);

        xTs = C[I[i]] * xTs;

        sparse_operator::axpy(xTs, xi, Hs, w_size);
    }
    for (i = 0; i < w_size; i++) Hs[i] = s[i] + 2 * Hs[i];
    if (regularize_bias == 0) Hs[w_size - 1] -= s[w_size - 1];
}

void l2r_l2_svc_fun::subXTv(double *v, double *XTv) {
    int w_size = get_nr_variable();
    std::fill(XTv, XTv + w_size, 0);
    for (int i = 0; i < sizeI; i++) {
        sparse_operator::axpy(v[i], prob->m_x[I[i]], XTv, w_size);
    }
}

class l2r_l2_svr_fun : public l2r_l2_svc_fun {
public:
    l2r_l2_svr_fun(const problem *prob, const parameter *param, double *C);

    void grad(double *w, double *g) override;

private:
    double C_times_loss(int i, double wx_i) override;
    double p;
};

l2r_l2_svr_fun::l2r_l2_svr_fun(const problem *prob, const parameter *param, double *C)
    : l2r_l2_svc_fun(prob, param, C) {
    this->p = param->p;
    this->regularize_bias = param->regularize_bias;
}

double l2r_l2_svr_fun::C_times_loss(int i, double wx_i) {
    double d = wx_i - prob->m_y[i];
    if (d < -p)
        return C[i] * (d + p) * (d + p);
    else if (d > p)
        return C[i] * (d - p) * (d - p);
    return 0;
}

void l2r_l2_svr_fun::grad(double *w, double *g) {
    int i;
    auto &y = prob->m_y;
    int l = prob->l;
    int w_size = get_nr_variable();
    double d;

    sizeI = 0;
    for (i = 0; i < l; i++) {
        d = wx[i] - y[i];

        // generate index set I
        if (d < -p) {
            tmp[sizeI] = C[i] * (d + p);
            I[sizeI] = i;
            sizeI++;
        } else if (d > p) {
            tmp[sizeI] = C[i] * (d - p);
            I[sizeI] = i;
            sizeI++;
        }
    }
    subXTv(tmp.data(), g);

    for (i = 0; i < w_size; i++) g[i] = w[i] + 2 * g[i];
    if (regularize_bias == 0) g[w_size - 1] -= w[w_size - 1];
}

static int solve_l2r_l1l2_svr(const problem *prob, const parameter *param, double *w, int max_iter = 300) {
    const int solver_type = param->solver_type;
    int l = prob->l;
    double C = param->C;
    double p = param->p;
    int w_size = prob->n;
    double eps = param->eps;
    int i, s, iter = 0;
    int active_size = l;
    std::vector<int> index(l, 0);

    double d, G, H;
    double Gmax_old = HUGE_VAL;
    double Gmax_new, Gnorm1_new;
    double Gnorm1_init = -1.0;  // Gnorm1_init is initialized at the first iteration
    std::vector<double> beta(l, NAN);
    std::vector<double> QD(l, NAN);
    auto &y = prob->m_y;

    // L2R_L2LOSS_SVR_DUAL
    double lambda = 0.5 / C;
    double upper_bound = HUGE_VAL;

    if (solver_type == L2R_L1LOSS_SVR_DUAL) {
        lambda = 0;
        upper_bound = C;
    }

    // Initial beta can be set here. Note that
    // -upper_bound <= beta[i] <= upper_bound
    for (i = 0; i < l; i++) beta[i] = 0;

    for (i = 0; i < w_size; i++) w[i] = 0;
    for (i = 0; i < l; i++) {
        auto xi = prob->m_x[i];
        QD[i] = sparse_operator::nrm2_sq(xi, w_size);
        sparse_operator::axpy(beta[i], xi, w, w_size);

        index[i] = i;
    }

    while (iter < max_iter) {
        Gmax_new = 0;
        Gnorm1_new = 0;

        for (i = 0; i < active_size; i++) {
            int j = i + rand() % (active_size - i);
            std::swap(index[i], index[j]);
        }

        for (s = 0; s < active_size; s++) {
            i = index[s];
            G = -y[i] + lambda * beta[i];
            H = QD[i] + lambda;

            auto xi = prob->m_x[i];
            G += sparse_operator::dot(w, xi, w_size);

            double Gp = G + p;
            double Gn = G - p;
            double violation = 0;
            if (beta[i] == 0) {
                if (Gp < 0)
                    violation = -Gp;
                else if (Gn > 0)
                    violation = Gn;
                else if (Gp > Gmax_old && Gn < -Gmax_old) {
                    active_size--;
                    std::swap(index[s], index[active_size]);
                    s--;
                    continue;
                }
            } else if (beta[i] >= upper_bound) {
                if (Gp > 0)
                    violation = Gp;
                else if (Gp < -Gmax_old) {
                    active_size--;
                    std::swap(index[s], index[active_size]);
                    s--;
                    continue;
                }
            } else if (beta[i] <= -upper_bound) {
                if (Gn < 0)
                    violation = -Gn;
                else if (Gn > Gmax_old) {
                    active_size--;
                    std::swap(index[s], index[active_size]);
                    s--;
                    continue;
                }
            } else if (beta[i] > 0)
                violation = fabs(Gp);
            else
                violation = fabs(Gn);

            Gmax_new = std::max(Gmax_new, violation);
            Gnorm1_new += violation;

            // obtain Newton direction d
            if (Gp < H * beta[i])
                d = -Gp / H;
            else if (Gn > H * beta[i])
                d = -Gn / H;
            else
                d = -beta[i];

            if (fabs(d) < 1.0e-12) continue;

            double beta_old = beta[i];
            beta[i] = std::min(std::max(beta[i] + d, -upper_bound), upper_bound);
            d = beta[i] - beta_old;

            if (d != 0) sparse_operator::axpy(d, xi, w, w_size);
        }

        if (iter == 0) Gnorm1_init = Gnorm1_new;
        iter++;
        if (iter % 10 == 0) printf(".");

        if (Gnorm1_new <= eps * Gnorm1_init) {
            if (active_size == l)
                break;
            else {
                active_size = l;
                printf("*");
                Gmax_old = HUGE_VAL;
                continue;
            }
        }

        Gmax_old = Gmax_new;
    }

    printf("\noptimization finished, #iter = %d\n", iter);

    // calculate objective value
    double v = 0;
    int nSV = 0;
    for (i = 0; i < w_size; i++) v += w[i] * w[i];
    v = 0.5 * v;
    for (i = 0; i < l; i++) {
        v += p * fabs(beta[i]) - y[i] * beta[i] + 0.5 * lambda * beta[i] * beta[i];
        if (beta[i] != 0) nSV++;
    }

    printf("Objective value = %lf\n", v);
    printf("nSV = %d\n", nSV);
    return iter;
}

static void train_one(const problem *prob, const parameter *param, double *w, double Cp, double Cn) {
    int solver_type = param->solver_type;
    int dual_solver_max_iter = 300;
    int iter;

    bool is_regression =
        (solver_type == L2R_L2LOSS_SVR || solver_type == L2R_L1LOSS_SVR_DUAL || solver_type == L2R_L2LOSS_SVR_DUAL);

    // Some solvers use Cp,Cn but not C array; extensions possible but no plan for now
    double *C = new double[prob->l];
    double primal_solver_tol = param->eps;
    if (is_regression) {
        for (int i = 0; i < prob->l; i++) C[i] = param->C;
    } else {
        int pos = 0;
        for (int i = 0; i < prob->l; i++) {
            if (prob->m_y[i] > 0) {
                pos++;
                C[i] = Cp;
            } else
                C[i] = Cn;
        }
        int neg = prob->l - pos;
        primal_solver_tol = param->eps * std::max(std::min(pos, neg), 1) / prob->l;
    }

    switch (solver_type) {
        case L2R_L2LOSS_SVR: {
            l2r_l2_svr_fun fun_obj(prob, param, C);
            NEWTON newton_obj(&fun_obj, primal_solver_tol);
            newton_obj.newton(w);
            break;
        }
        case L2R_L1LOSS_SVR_DUAL: {
            iter = solve_l2r_l1l2_svr(prob, param, w, dual_solver_max_iter);
            if (iter >= dual_solver_max_iter)
                printf("\nWARNING: reaching max number of iterations\nUsing -s 11 may be faster (also see FAQ)\n\n");

            break;
        }
        case L2R_L2LOSS_SVR_DUAL: {
            iter = solve_l2r_l1l2_svr(prob, param, w, dual_solver_max_iter);
            if (iter >= dual_solver_max_iter) {
                printf("\nWARNING: reaching max number of iterations\nSwitching to use -s 11\n\n");
                // primal_solver_tol obtained from eps for dual may be too loose
                primal_solver_tol *= 0.001;
                l2r_l2_svr_fun fun_obj(prob, param, C);
                NEWTON newton_obj(&fun_obj, primal_solver_tol);
                newton_obj.newton(w);
            }
            break;
        }
        default:
            fprintf(stderr, "ERROR: unknown solver_type\n");
            break;
    }

    delete[] C;
}

// Calculate the initial C for parameter selection
static double calc_start_C(const problem *prob, const parameter *param) {
    int i;
    double xTx, max_xTx;
    max_xTx = 0;
    for (i = 0; i < prob->l; i++) {
        xTx = 0;
        auto *xi = prob->m_x[i];
        for (int j = 0; j < prob->n; ++j) {
            xTx += xi[j] * xi[j];
        }
        if (xTx > max_xTx) max_xTx = xTx;
    }

    double min_C = 1.0;
    if (param->solver_type == L2R_L2LOSS_SVR) {
        double sum_y, loss, y_abs;
        double delta2 = 0.1;
        sum_y = 0, loss = 0;
        for (i = 0; i < prob->l; i++) {
            y_abs = fabs(prob->m_y[i]);
            sum_y += y_abs;
            loss += std::max(y_abs - param->p, 0.0) * std::max(y_abs - param->p, 0.0);
        }
        if (loss > 0)
            min_C = delta2 * delta2 * loss / (8 * sum_y * sum_y * max_xTx);
        else
            min_C = HUGE_VAL;
    }

    return pow(2, floor(log(min_C) / log(2.0)));
}

static double calc_max_p(const problem *prob) {
    int i;
    double max_p = 0.0;
    for (i = 0; i < prob->l; i++) max_p = std::max(max_p, fabs(prob->m_y[i]));

    return max_p;
}

static void find_parameter_C(const problem *prob, parameter *param_tmp, double start_C, double max_C, double *best_C,
                             double *best_score, const int *fold_start, const int *perm, const problem *subprob,
                             int nr_fold) {
    // variables for CV
    std::vector<double> target(prob->l);

    // variables for warm start
    double ratio = 2;
    std::vector<std::vector<double>> prev_w(nr_fold);
    int num_unchanged_w = 0;

    if (param_tmp->solver_type == L2R_L2LOSS_SVR) *best_score = HUGE_VAL;
    *best_C = start_C;

    param_tmp->C = start_C;
    while (param_tmp->C <= max_C) {
        for (int i = 0; i < nr_fold; i++) {
            int j;
            int begin = fold_start[i];
            int end = fold_start[i + 1];

            param_tmp->init_sol = prev_w[i];
            model submodel = train(&subprob[i], param_tmp);

            int total_w_size;
            total_w_size = subprob[i].n;

            if (prev_w[i].empty()) {
                prev_w[i] = submodel.w;
            } else if (num_unchanged_w >= 0) {
                double norm_w_diff = 0;
                for (j = 0; j < total_w_size; j++) {
                    norm_w_diff += (submodel.w[j] - prev_w[i][j]) * (submodel.w[j] - prev_w[i][j]);
                    prev_w[i][j] = submodel.w[j];
                }
                norm_w_diff = sqrt(norm_w_diff);

                if (norm_w_diff > 1e-15) num_unchanged_w = -1;
            } else {
                prev_w[i] = submodel.w;
            }

            for (j = begin; j < end; j++) target[perm[j]] = predict(&submodel, prob->m_x[perm[j]]);
        }

        if (param_tmp->solver_type == L2R_L2LOSS_SVR) {
            double total_error = 0.0;
            for (int i = 0; i < prob->l; i++) {
                double y = prob->m_y[i];
                double v = target[i];
                total_error += (v - y) * (v - y);
            }
            double current_error = total_error / prob->l;
            if (current_error < *best_score) {
                *best_C = param_tmp->C;
                *best_score = current_error;
            }

            if (param_tmp->verbose) {
                printf("log2c=%7.2f\tp=%7.2f\tMean squared error=%g\n", log(param_tmp->C) / log(2.0),
                       param_tmp->p, current_error);
            }
        }

        num_unchanged_w++;
        if (num_unchanged_w == 5) break;
        param_tmp->C = param_tmp->C * ratio;
    }

    if (param_tmp->C > max_C) printf("WARNING: maximum C reached.\n");
}

//
// Interface functions
//
model train(const problem *prob, const parameter *param) {
    int n = prob->n;
    int w_size = prob->n;
    model model_;

    if (prob->bias >= 0)
        model_.nr_feature = n - 1;
    else
        model_.nr_feature = n;
    model_.param = *param;
    model_.bias = prob->bias;

    if (check_regression_model(&model_)) {
        model_.w.resize(w_size);

        if (!param->init_sol.empty())
            model_.w = param->init_sol;
        else
            std::fill(model_.w.begin(), model_.w.end(), 0);

        train_one(prob, param, model_.w.data(), 0, 0);
    }
    return model_;
}

void cross_validation(const problem *prob, const parameter *param, int nr_fold, double *target) {
    int l = prob->l;
    std::vector<int> perm(l);
    if (nr_fold > l) {
        nr_fold = l;
        fprintf(
            stderr,
            "WARNING: # folds > # data. Will use # folds = # data instead (i.e., leave-one-out cross validation)\n");
    }
    std::vector<int> fold_start(nr_fold + 1);
    for (int i = 0; i < l; i++) perm[i] = i;
    for (int i = 0; i < l; i++) {
        int j = i + rand() % (l - i);
        std::swap(perm[i], perm[j]);
    }
    for (int i = 0; i <= nr_fold; i++) fold_start[i] = i * l / nr_fold;

    for (int i = 0; i < nr_fold; i++) {
        int begin = fold_start[i];
        int end = fold_start[i + 1];
        int j, k;
        struct problem subprob;

        subprob.bias = prob->bias;
        subprob.n = prob->n;
        subprob.l = l - (end - begin);
        subprob.m_x.resize(subprob.l);
        subprob.m_y.resize(subprob.l);

        k = 0;
        for (j = 0; j < begin; j++) {
            subprob.m_x[k] = prob->m_x[perm[j]];
            subprob.m_y[k] = prob->m_y[perm[j]];
            ++k;
        }
        for (j = end; j < l; j++) {
            subprob.m_x[k] = prob->m_x[perm[j]];
            subprob.m_y[k] = prob->m_y[perm[j]];
            ++k;
        }
        model submodel = train(&subprob, param);
        for (j = begin; j < end; j++) target[perm[j]] = predict(&submodel, prob->m_x[perm[j]]);
    }
}

void find_parameters(const problem *prob, const parameter *param, int nr_fold, double start_C, double start_p,
                     double *best_C, double *best_p, double *best_score) {
    int l = prob->l;
    std::vector<int> perm(l);
    std::vector<problem> subprob(nr_fold);

    if (nr_fold > l) {
        nr_fold = l;
        fprintf(
            stderr,
            "WARNING: # folds > # data. Will use # folds = # data instead (i.e., leave-one-out cross validation)\n");
    }
    std::vector<int> fold_start(nr_fold + 1);
    for (int i = 0; i < l; i++) perm[i] = i;
    for (int i = 0; i < l; i++) {
        int j = i + rand() % (l - i);
        std::swap(perm[i], perm[j]);
    }
    for (int i = 0; i <= nr_fold; i++) fold_start[i] = i * l / nr_fold;

    for (int i = 0; i < nr_fold; i++) {
        int begin = fold_start[i];
        int end = fold_start[i + 1];
        int j, k;

        subprob[i].bias = prob->bias;
        subprob[i].n = prob->n;
        subprob[i].l = l - (end - begin);
        subprob[i].m_x.resize(subprob[i].l);
        subprob[i].m_y.resize(subprob[i].l);

        k = 0;
        for (j = 0; j < begin; j++) {
            subprob[i].m_x[k] = prob->m_x[perm[j]];
            subprob[i].m_y[k] = prob->m_y[perm[j]];
            ++k;
        }
        for (j = end; j < l; j++) {
            subprob[i].m_x[k] = prob->m_x[perm[j]];
            subprob[i].m_y[k] = prob->m_y[perm[j]];
            ++k;
        }
    }

    struct parameter param_tmp = *param;
    *best_p = -1;
    if (param->solver_type == L2R_L2LOSS_SVR) {
        double max_p = calc_max_p(prob);
        int num_p_steps = 20;
        double max_C = 1048576;
        *best_score = HUGE_VAL;

        int i = num_p_steps - 1;
        if (start_p > 0) i = std::min((int)(start_p / (max_p / num_p_steps)), i);
        for (; i >= 0; i--) {
            param_tmp.p = i * max_p / num_p_steps;
            double start_C_tmp;
            if (start_C <= 0)
                start_C_tmp = calc_start_C(prob, &param_tmp);
            else
                start_C_tmp = start_C;
            start_C_tmp = std::min(start_C_tmp, max_C);
            double best_C_tmp, best_score_tmp;

            find_parameter_C(prob, &param_tmp, start_C_tmp, max_C, &best_C_tmp, &best_score_tmp, fold_start.data(),
                             perm.data(), subprob.data(), nr_fold);

            if (best_score_tmp < *best_score) {
                *best_p = param_tmp.p;
                *best_C = best_C_tmp;
                *best_score = best_score_tmp;
            }
        }
    }
}

double predict_values(const struct model *model_, const double *x) {
    int n = model_->nr_feature;
    if (model_->bias >= 0) n = model_->nr_feature + 1;
    const double *w = model_->w.data();

    double dec_value = 0;
    for (int i = 0; i < n; ++i) {
        dec_value += w[i] * x[i];
    }
    return dec_value;
}

double predict(const model *model_, const double *x) { return predict_values(model_, x); }

static const char *solver_type_table[] = {
    "", "", "", "", "", "", "", "", "",     "", "", "L2R_L2LOSS_SVR", "L2R_L2LOSS_SVR_DUAL", "L2R_L1LOSS_SVR_DUAL",
    "", "", "", "", "", "", "", "", nullptr};

int save_model(const char *model_file_name, const struct model *model_, const std::vector<std::string>& f_names) {
    int nr_feature = model_->nr_feature;
    const parameter &param = model_->param;

    FILE *fp = fopen(model_file_name, "w");
    if (fp == nullptr) return -1;

    fprintf(fp, "solver_type,%s\n", solver_type_table[param.solver_type]);
    fprintf(fp, "nr_feature,%d\n", nr_feature);
    fprintf(fp, "bias,%.17g\n", model_->bias);
    for (int i = 0; i < nr_feature; i++) {
        fprintf(fp, "%s,", f_names[i].c_str());
        fprintf(fp, "%.17g", model_->w[i]);
        fprintf(fp, "\n");
    }
    if (model_->bias >= 0) {
        fprintf(fp, "INTERCEPT,%.17g\n", model_->w[nr_feature]);
    }

    if (ferror(fp) != 0 || fclose(fp) != 0)
        return -1;
    else
        return 0;
}

model *load_model(const char *model_file_name, std::vector<std::string>& f_names) {
    f_names.clear();
    std::ifstream ifs(model_file_name, std::ifstream::in);

    if (!ifs) {
        printf("load model open file %s failed\n", model_file_name);
        return nullptr;
    }

    auto *model_ = new model;
    parameter &param = model_->param;

    std::string s;
    while (getline(ifs, s)) {
        if (s.empty() || s.front() == '#') continue;
        auto itr = s.find_first_of(',');
        auto first = s.substr(0, itr);
        auto second = s.substr(itr + 1);
        if (first == "nr_feature") {
            model_->nr_feature = std::stoi(second);
        } else if (first == "bias") {
            model_->bias = std::stod(second);
        } else if (first == "solver_type") {
            for (int i = 0; solver_type_table[i]; i++) {
                if (strcmp(solver_type_table[i], second.c_str()) == 0) {
                    param.solver_type = i;
                    break;
                }
            }
        } else if (first == "INTERCEPT") {
            model_->w.push_back(std::stod(second));
        } else {
            f_names.push_back(first);
            model_->w.push_back(std::stod(second));
        }
    }

    ifs.close();
    return model_;
}

const char *check_parameter(const problem *prob, const parameter *param) {
    if (param->eps <= 0) return "eps <= 0";

    if (param->C <= 0) return "C <= 0";

    if (param->p < 0) return "p < 0";

    if (!param->regularize_bias) {
        if (prob->bias != 1.0) return "To not regularize bias, must specify -B 1 along with -R";
        if (param->solver_type != L2R_L2LOSS_SVR) return "-R option supported only for solver L2R_L2LOSS_SVR";
    }

    if (param->solver_type != L2R_L2LOSS_SVR && param->solver_type != L2R_L2LOSS_SVR_DUAL &&
        param->solver_type != L2R_L1LOSS_SVR_DUAL)
        return "unknown solver type";

    if (!param->init_sol.empty() && param->solver_type != L2R_L2LOSS_SVR)
        return "Initial-solution specification supported only for solvers L2R_LR, L2R_L2LOSS_SVC";

    return nullptr;
}

int check_regression_model(const struct model *model_) {
    return (model_->param.solver_type == L2R_L2LOSS_SVR || model_->param.solver_type == L2R_L1LOSS_SVR_DUAL ||
            model_->param.solver_type == L2R_L2LOSS_SVR_DUAL);
}

// On entry *f must be the function value of w
// On exit w is updated and *f is the new function value
double function::linesearch_and_update(double *w, double *s, double *f, double *g, double alpha) {
    double gTs = 0;
    double eta = 0.01;
    int n = get_nr_variable();
    int max_num_linesearch = 20;
    double *w_new = new double[n];
    double fold = *f;

    for (int i = 0; i < n; i++) gTs += s[i] * g[i];

    int num_linesearch = 0;
    for (num_linesearch = 0; num_linesearch < max_num_linesearch; num_linesearch++) {
        for (int i = 0; i < n; i++) w_new[i] = w[i] + alpha * s[i];
        *f = fun(w_new);
        if (*f - fold <= eta * alpha * gTs)
            break;
        else
            alpha *= 0.5;
    }

    if (num_linesearch >= max_num_linesearch) {
        *f = fold;
        return 0;
    } else
        memcpy(w, w_new, sizeof(double) * n);

    delete[] w_new;
    return alpha;
}

NEWTON::NEWTON(const function *fun_obj, double eps, double eps_cg, int max_iter) {
    this->fun_obj = const_cast<function *>(fun_obj);
    this->eps = eps;
    this->eps_cg = eps_cg;
    this->max_iter = max_iter;
}

void NEWTON::newton(double *w) {
    int n = fun_obj->get_nr_variable();
    int i;
    double step_size;
    double f, fold, actred;
    double init_step_size = 1;
    int search = 1, iter = 1, inc = 1;
    double *s = new double[n];
    double *r = new double[n];
    double *g = new double[n];

    const double alpha_pcg = 0.01;
    double *M = new double[n];

    // calculate gradient norm at w=0 for stopping condition.
    double *w0 = new double[n];
    for (i = 0; i < n; i++) w0[i] = 0;
    fun_obj->fun(w0);
    fun_obj->grad(w0, g);
    double gnorm0 = dnrm2_(&n, g, &inc);
    delete[] w0;

    f = fun_obj->fun(w);
    fun_obj->grad(w, g);
    double gnorm = dnrm2_(&n, g, &inc);
    if (verbose) printf("init f %5.3e |g| %5.3e\n", f, gnorm);

    if (gnorm <= eps * gnorm0) search = 0;

    while (iter <= max_iter && search) {
        fun_obj->get_diag_preconditioner(M);
        for (i = 0; i < n; i++) M[i] = (1 - alpha_pcg) + alpha_pcg * M[i];
        pcg(g, M, s, r);

        fold = f;
        step_size = fun_obj->linesearch_and_update(w, s, &f, g, init_step_size);

        if (step_size == 0) {
            if (verbose) printf("WARNING: line search fails\n");
            break;
        }

        fun_obj->grad(w, g);
        gnorm = dnrm2_(&n, g, &inc);

        // printf("iter %2d f %5.3e |g| %5.3e CG %3d step_size %4.2e \n", iter, f, gnorm, cg_iter, step_size);

        if (gnorm <= eps * gnorm0) break;
        if (f < -1.0e+32) {
            if (verbose) printf("WARNING: f < -1.0e+32\n");
            break;
        }
        actred = fold - f;
        if (fabs(actred) <= 1.0e-12 * fabs(f)) {
            if (verbose) printf("WARNING: actred too small\n");
            break;
        }

        iter++;
    }

    if (iter >= max_iter) printf("\nWARNING: reaching max number of Newton iterations\n");

    delete[] g;
    delete[] r;
    delete[] s;
    delete[] M;
}

int NEWTON::pcg(double *g, double *M, double *s, double *r) {
    int i, inc = 1;
    int n = fun_obj->get_nr_variable();
    double one = 1;
    double *d = new double[n];
    double *Hd = new double[n];
    double zTr, znewTrnew, alpha, beta, cgtol, dHd;
    double *z = new double[n];
    double Q = 0, newQ, Qdiff;

    for (i = 0; i < n; i++) {
        s[i] = 0;
        r[i] = -g[i];
        z[i] = r[i] / M[i];
        d[i] = z[i];
    }

    zTr = ddot_(&n, z, &inc, r, &inc);
    double gMinv_norm = sqrt(zTr);
    cgtol = std::min(eps_cg, sqrt(gMinv_norm));
    int cg_iter = 0;
    int max_cg_iter = std::max(n, 5);

    while (cg_iter < max_cg_iter) {
        cg_iter++;

        fun_obj->Hv(d, Hd);
        dHd = ddot_(&n, d, &inc, Hd, &inc);
        // avoid 0/0 in getting alpha
        if (dHd <= 1.0e-16) break;

        alpha = zTr / dHd;
        daxpy_(&n, &alpha, d, &inc, s, &inc);
        alpha = -alpha;
        daxpy_(&n, &alpha, Hd, &inc, r, &inc);

        // Using quadratic approximation as CG stopping criterion
        newQ = -0.5 * (ddot_(&n, s, &inc, r, &inc) - ddot_(&n, s, &inc, g, &inc));
        Qdiff = newQ - Q;
        if (newQ <= 0 && Qdiff <= 0) {
            if (cg_iter * Qdiff >= cgtol * newQ) break;
        } else {
            printf("WARNING: quadratic approximation > 0 or increasing in CG\n");
            break;
        }
        Q = newQ;

        for (i = 0; i < n; i++) z[i] = r[i] / M[i];
        znewTrnew = ddot_(&n, z, &inc, r, &inc);
        beta = znewTrnew / zTr;
        dscal_(&n, &beta, d, &inc);
        daxpy_(&n, &one, z, &inc, d, &inc);
        zTr = znewTrnew;
    }

    if (cg_iter == max_cg_iter) printf("WARNING: reaching maximal number of CG steps\n");

    delete[] d;
    delete[] Hd;
    delete[] z;

    return cg_iter;
}

double LinearSvm::predict_x(const double* x) {
    return predict(&model_, x);
}

void LinearSvm::predict_test() {
    std::vector<double> predicts(prob.l);
    for (int i = 0; i < prob.l; ++i) {
        predicts[i] = predict(&model_, prob.m_x[i]);
    }

    int correct = 0;
    int total = 0;
    double error = 0;
    double sump = 0, sumt = 0, sumpp = 0, sumtt = 0, sumpt = 0;
    for (int i = 0; i < prob.l; ++i) {
        double target_label = prob.m_y[i];
        double predict_label = predicts[i];

        if (predict_label == target_label) ++correct;
        error += (predict_label - target_label) * (predict_label - target_label);
        sump += predict_label;
        sumt += target_label;
        sumpp += predict_label * predict_label;
        sumtt += target_label * target_label;
        sumpt += predict_label * target_label;
        ++total;
    }

    printf("Mean squared error = %g (regression)\n", error / total);
    printf("R2 = %g (regression)\n",
           ((total * sumpt - sump * sumt) * (total * sumpt - sump * sumt)) /
               ((total * sumpp - sump * sump) * (total * sumtt - sumt * sumt)));
}

bool LinearSvm::load(std::string model_file_name, std::vector<std::string>& f_names) {
    model* model_tmp = load_model(model_file_name.c_str(), f_names);

    if (model_tmp) {
        model_ = *model_tmp;
        delete model_tmp;
        return true;
    } else {
        return false;
    }
}

bool LinearSvm::save(std::string model_file_name, const std::vector<std::string>& f_names) {
    if (save_model(model_file_name.c_str(), &model_, f_names)) {
        fprintf(stderr, "can't save model to file %s\n", model_file_name.c_str());
        return false;
    } else {
        printf("save model to %s\n", model_file_name.c_str());
        return true;
    }
}

void LinearSvm::train_model_with_param() {
    model_ = train(&prob, &param);
}

void LinearSvm::work() {
    const char *error_msg = check_parameter(&prob, &param);

    if (error_msg) {
        fprintf(stderr, "ERROR: %s\n", error_msg);
        exit(1);
    }

    if (flag_find_parameters) {
        do_find_parameters();
        train_model_with_param();
    } else if (flag_cross_validation) {
        do_cross_validation();
    } else {
        train_model_with_param();
    }
}

void LinearSvm::set_data(ornate::DataSet &ds) {
    prob.m_x.clear();
    prob.m_y.swap(ds.m_y);
    for (size_t i = 0; i < prob.m_y.size(); ++i) {
        prob.m_x.push_back(ds.m_xs.data() + ds.m_f_count * i);
    }
    prob.n = (int)ds.m_f_count;
    prob.l = (int)prob.m_y.size();
    prob.bias = bias;
}

void LinearSvm::do_find_parameters() {
    double start_C, start_p, best_C, best_p, best_score;
    if (flag_C_specified)
        start_C = param.C;
    else
        start_C = -1.0;
    if (flag_p_specified)
        start_p = param.p;
    else
        start_p = -1.0;

    printf("Doing parameter search with %d-fold cross validation.\n", nr_fold);
    find_parameters(&prob, &param, nr_fold, start_C, start_p, &best_C, &best_p, &best_score);
    if (param.solver_type == L2R_L2LOSS_SVR)
        printf("Best C = %g Best p = %g  CV MSE = %g\n", best_C, best_p, best_score);

    param.C = best_C;
    param.p = best_p;
}

void LinearSvm::do_cross_validation() {
    int i;
    double total_error = 0;
    double sumv = 0, sumy = 0, sumvv = 0, sumyy = 0, sumvy = 0;
    std::vector<double> target(prob.l);

    cross_validation(&prob, &param, nr_fold, target.data());
    if (param.solver_type == L2R_L2LOSS_SVR || param.solver_type == L2R_L1LOSS_SVR_DUAL ||
        param.solver_type == L2R_L2LOSS_SVR_DUAL) {
        for (i = 0; i < prob.l; i++) {
            double y = prob.m_y[i];
            double v = target[i];
            total_error += (v - y) * (v - y);
            sumv += v;
            sumy += y;
            sumvv += v * v;
            sumyy += y * y;
            sumvy += v * y;
        }
        printf("Cross Validation Mean squared error = %g\n", total_error / prob.l);
        printf("Cross Validation Squared correlation coefficient = %g\n",
               ((prob.l * sumvy - sumv * sumy) * (prob.l * sumvy - sumv * sumy)) /
                   ((prob.l * sumvv - sumv * sumv) * (prob.l * sumyy - sumy * sumy)));
    }
}

void LinearSvm::init() {
    // default solver for parameter selection is L2R_L2LOSS_SVC
    if (flag_find_parameters) {
        if (!flag_cross_validation) nr_fold = 5;
        if (!flag_solver_specified) {
            fprintf(stderr, "Solver not specified. Using -s 11\n");
            param.solver_type = L2R_L2LOSS_SVR;
        } else if (param.solver_type != L2R_L2LOSS_SVR) {
            fprintf(stderr, "Warm-start parameter search only available for -s 0, -s 2 and -s 11\n");
            exit_with_help();
        }
    }

    if (param.eps == HUGE_VAL) {
        switch (param.solver_type) {
            case L2R_L2LOSS_SVR:
                param.eps = 0.0001;
                break;
            case L2R_L1LOSS_SVR_DUAL:
            case L2R_L2LOSS_SVR_DUAL:
                param.eps = 0.1;
                break;
        }
    }
}

void LinearSvm::exit_with_help() {
    printf(
        "Usage: train [options] training_set_file [model_file]\n"
        "options:\n"
        "-s type : set type of solver (default 1)\n"
        "  for regression\n"
        "       11 -- L2-regularized L2-loss support vector regression (primal)\n"
        "       12 -- L2-regularized L2-loss support vector regression (dual)\n"
        "       13 -- L2-regularized L1-loss support vector regression (dual)\n"
        "-c cost : set the parameter C (default 1)\n"
        "-p epsilon : set the epsilon in loss function of SVR (default 0.1)\n"
        "-n nu : set the parameter nu of one-class SVM (default 0.5)\n"
        "-e epsilon : set tolerance of termination criterion\n"
        "       -s 11\n"
        "               |f'(w)|_2 <= eps*|f'(w0)|_2 (default 0.0001)\n"
        "       -s 12 and 13\n"
        "               |f'(alpha)|_1 <= eps |f'(alpha0)|,\n"
        "               where f is the dual function (default 0.1)\n"
        "-B bias : if bias >= 0, instance x becomes [x; bias]; if < 0, no bias term added (default -1)\n"
        "-R : not regularize the bias; must with -B 1 to have the bias; DON'T use this unless you know what it is\n"
        "       (for -s 0, 2, 5, 6, 11)\n"
        "-v n: n-fold cross validation mode\n"
        "-C : find parameters (C for -s 0, 2 and C, p for -s 11)\n");
    exit(1);
}

}  // namespace ornate::svm
