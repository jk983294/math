#ifndef _LIBLINEAR_H
#define _LIBLINEAR_H

#include <math_ds.h>
#include <cmath>
#include <string>
#include <vector>

namespace ornate::svm {

struct problem {
    int l;  // the number of training data
    int n;  // the number of feature (including the bias feature if bias >= 0)
    std::vector<double> m_y;
    std::vector<const double *> m_x;
    double bias{1}; /* < 0 if no bias term */
};

enum { L2R_L2LOSS_SVR = 11, L2R_L2LOSS_SVR_DUAL, L2R_L1LOSS_SVR_DUAL }; /* solver_type */

struct parameter {
    int solver_type = L2R_L2LOSS_SVR;

    /* these are for training only */
    double eps{HUGE_VAL};          /* stopping tolerance */
    double C{1};                   // cost of constraints violation
    double p{0.1};                 // sensitiveness of loss of support vector regression
    std::vector<double> init_sol;  // initial weight vectors
    bool regularize_bias{true};
};

struct model {
    struct parameter param;
    int nr_feature{0}; /* number of features */
    std::vector<double> w;
    double bias;
};

class function {
public:
    virtual double fun(double *w) = 0;
    virtual void grad(double *w, double *g) = 0;
    virtual void Hv(double *s, double *Hs) = 0;
    virtual int get_nr_variable() = 0;
    virtual void get_diag_preconditioner(double *M) = 0;
    virtual ~function() = default;

    // base implementation in newton.cpp
    virtual double linesearch_and_update(double *w, double *s, double *f, double *g, double alpha);
};

class NEWTON {
public:
    NEWTON(const function *fun_obj, double eps = 0.1, double eps_cg = 0.5, int max_iter = 1000);
    ~NEWTON() = default;

    void newton(double *w);

private:
    int pcg(double *g, double *M, double *s, double *r);

    double eps;
    double eps_cg;
    int max_iter;
    function *fun_obj;
};

struct model train(const struct problem *prob, const struct parameter *param);
void cross_validation(const struct problem *prob, const struct parameter *param, int nr_fold, double *target);
void find_parameters(const struct problem *prob, const struct parameter *param, int nr_fold, double start_C,
                     double start_p, double *best_C, double *best_p, double *best_score);

double predict_values(const struct model *model_, const double *x);
double predict(const struct model *model_, const double *x);

int save_model(const char *model_file_name, const struct model *model_);
model *load_model(const char *model_file_name);

const char *check_parameter(const struct problem *prob, const struct parameter *param);
int check_regression_model(const struct model *model);

struct LinearSvm {
    static void exit_with_help();
    void init();
    void work();
    void do_cross_validation();
    void do_find_parameters();
    void set_data(ornate::DataSet &ds);
    void predict_test();

private:
    void train_model_with_param();

public:
    struct parameter param;
    struct problem prob;
    struct model model_;
    bool flag_cross_validation{false};
    bool flag_find_parameters{false};
    bool flag_C_specified{false};
    bool flag_p_specified{false};
    bool flag_solver_specified{false};
    int nr_fold{0};
    double bias{1};
    std::string model_file_name;
};

}  // namespace ornate::svm

#endif /* _LIBLINEAR_H */
