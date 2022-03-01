#include <getopt.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include "linear_svm.h"

using namespace ornate::svm;

void prepare_data(ornate::DataSet& ds);

int main(int argc, char** argv) {
    LinearSvm lsvm;
    int opt;
    while ((opt = getopt(argc, argv, "hs:c:p:e:B:v:CR")) != -1) {
        switch (opt) {
            case 's':
                lsvm.param.solver_type = std::stoi(optarg);
                break;
            case 'c':
                lsvm.param.C = std::stod(optarg);
                lsvm.flag_C_specified = true;
                break;
            case 'p':
                lsvm.param.p = std::stod(optarg);
                lsvm.flag_p_specified = true;
                break;
            case 'e':
                lsvm.param.eps = std::stod(optarg);
                break;
            case 'B':
                lsvm.bias = std::stod(optarg);
                break;
            case 'v':
                lsvm.flag_cross_validation = true;
                lsvm.nr_fold = std::stoi(optarg);
                if (lsvm.nr_fold < 2) {
                    fprintf(stderr, "n-fold cross validation: n must >= 2\n");
                    lsvm.exit_with_help();
                }
                break;
            case 'C':
                lsvm.flag_find_parameters = true;
                break;
            case 'R':
                lsvm.param.regularize_bias = false;
                break;
            case 'h':
                lsvm.exit_with_help();
                return 0;
            default:
                abort();
        }
    }

    lsvm.init();
    ornate::DataSet ds;
    prepare_data(ds);
    lsvm.set_data(ds);
    lsvm.work();
    if (!lsvm.flag_cross_validation)
        lsvm.predict_test();
    return 0;
}

void prepare_data(ornate::DataSet& ds) {
    std::random_device rd;
    std::mt19937 generator(rd());
    std::normal_distribution<double> nd(0., 0.001);

    // y = 2 * x0 + 3 * x1  - 1 + delta
    size_t n = 100;
    ds.reserve(n * n, 3);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            double x0 = i, x1 = j, bias = -1;
            double y = 2 * x0 + 3 * x1 + bias + nd(generator);
            std::vector<double> xs = {x0, x1, 1};
            ds.add_row(y, xs);
        }
    }
}


