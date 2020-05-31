#ifndef ORNATE_L1TRENDFILTER_H
#define ORNATE_L1TRENDFILTER_H

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <numeric>
#include <vector>
#include "blas.h"
#include "lapack.h"

namespace ornate {

namespace {
const double done = 1.0;
const double dtwo = 2.0;
const double dminusone = -1.0;
const int ione = 1;
const int itwo = 2;
const int ithree = 3;
const int iseven = 7;
const double ALPHA = 0.01;    // linesearch parameter (0,0.5]
const double BETA = 0.5;      // linesearch parameter (0,1)
const double MU = 2;          // IPM parameter: t update
const double MAXITER = 40;    // IPM parameter: max iter. of IPM
const double MAXLSITER = 20;  // IPM parameter: max iter. of linesearch
const double TOL = 1e-4;      // IPM parameter: tolerance

/**
 * Computes y = D*x, where x has length n
 *     | 1 -2  1  0  0 |
 * y = | 0  1 -2  1  0 |*x
 *     | 0  0  1 -2  1 |
 * output y has n - 3 length
 */
void Dx(int n, const double *x, double *y) {
    for (int i = 0; i < n - 2; i++) y[i] = x[i] - x[i + 1] - x[i + 1] + x[i + 2];
}

/**
 * Computes y = D^T*x, where x has length n
 *     | 1  0  0 |
 *     |-2  1  0 |
 * y = | 1 -2  1 |*x
 *     | 0  1 -2 |
 *     | 0  0  1 |
 * output y has n + 2 length
 */
void DTx(int n, const double *x, double *y) {
    y[0] = x[0];
    y[1] = -x[0] - x[0] + x[1];
    for (int i = 2; i < n; i++) y[i] = x[i - 2] - x[i - 1] - x[i - 1] + x[i];  // y[2..n-1]
    y[n] = x[n - 2] - x[n - 1] - x[n - 1];                                     // y[n]
    y[n + 1] = x[n - 1];                                                       // y[n+1]
}

void yainvx(int n, const double a, const double *x, double *y) {
    while (n-- != 0) *y++ = a / *x++;
}
}  // namespace

struct L1TrendFilter {
    double l1tf_lambdamax(int n, double *y) {
        int m = n - 2;
        std::vector<double> vec(m, 0);
        std::vector<double> mat(7 * m, 0);
        std::vector<int> piv(m, 0);

        Dx(n, y, vec.data());
        double *dptr = mat.data();
        for (int i = 0; i < m; i++) {
            *dptr++ = 6;
            *dptr++ = -4.0;
            *dptr++ = 1.0;
        }

        /**
         * DPBSV computes the solution to system of linear equations A * X = B
         * The Cholesky decomposition is used to factor A
         * Lower triangle
         */
        int info;
        F77_CALL(dpbsv)("L", &m, &itwo, &ione, mat.data(), &ithree, vec.data(), &m, &info);
        if (info > 0) {  // if Cholesky factorization fails, try LU factorization
            dptr = mat.data();
            for (int i = 0; i < m; i++) {
                dptr++;
                dptr++;
                *dptr++ = 1.0;
                *dptr++ = -4.0;
                *dptr++ = 6.0;
                *dptr++ = -4.0;
                *dptr++ = 1.0;
            }

            /**
             * DGBSV computes the solution to system of linear equations A * X = B for GB matrices
             * The LU decomposition with partial pivoting and row interchanges is used
             */
            F77_CALL(dgbsv)(&m, &itwo, &itwo, &ione, mat.data(), &iseven, piv.data(), vec.data(), &m, &info);
            if (info > 0) return -1.0; /* if LU fails, return -1 */
        }
        double maxval = 0;
        for (int i = 0; i < m; i++) {
            if (fabs(vec[i]) > maxval) maxval = fabs(vec[i]);
        }
        return maxval;
    }

    void l1tf(const int n, const double *y, double *x) {
        const int m = n - 2;  // dimension length of Dx

        std::vector<double> S(3 * m, 0);  // matrix of size 3xm
        // DDTF stores Cholesky factor of DDT in packed symm.-band representation
        std::vector<double> DDTF(7 * m, 0);

        std::vector<double> DTz(n, 0);  // vector of size n
        std::vector<double> Dy(m, 0);   // vector of size m
        std::vector<double> DDTz(m, 0);
        std::vector<double> z(m, 0);      // dual variable
        std::vector<double> mu1(m, 1.0);  // dual of dual variables
        std::vector<double> mu2(m, 1.0);
        std::vector<double> f1(m, -lambda);  // constraints
        std::vector<double> f2(m, -lambda);
        std::vector<double> dz(m, 0);  // search directions
        std::vector<double> dmu1(m, 0);
        std::vector<double> dmu2(m, 0);
        std::vector<double> w(m, 0);
        std::vector<double> rz(m, 0);
        std::vector<double> tmp_m1(m, 0);
        std::vector<double> tmp_m2(m, 0);
        std::vector<int> IPIV(m, 0);

        /**
         * DDT : packed representation of symmetric banded matrix
         * fortran style storage (should be transposed in C)
         *  6  6  6  ...  6  6  6
         * -4 -4 -4  ... -4 -4  *
         *  1  1  1  ...  1  *  *
         */
        bool ddtf_chol = true;
        double *dptr = DDTF.data();
        for (int i = 0; i < m; i++) {
            *dptr++ = 6.0;
            *dptr++ = -4.0;
            *dptr++ = 1.0;
        }
        int info;
        F77_CALL(dpbtrf)("L", &m, &itwo, DDTF.data(), &ithree, &info);
        if (info > 0) {  // if Cholesky factorization fails, try LU factorization
            ddtf_chol = false;
            dptr = DDTF.data();
            for (int i = 0; i < m; i++) {
                dptr++;
                dptr++;
                *dptr++ = 1.0;
                *dptr++ = -4.0;
                *dptr++ = 6.0;
                *dptr++ = -4.0;
                *dptr++ = 1.0;
            }
            F77_CALL(dgbtrf)(&m, &m, &itwo, &itwo, DDTF.data(), &iseven, IPIV.data(), &info);
        }

        Dx(n, y, Dy.data());

        // barrier update parameter
        double t = -1, step = 1;

        DTx(m, z.data(), DTz.data());    // DTz = D'*z
        Dx(n, DTz.data(), DDTz.data());  // DDTz = D*D'*z

        if (is_debug) printf("%s %13s %12s %8s\n", "Iteration", "Primal obj.", "Dual obj.", "Gap");

        for (int iters = 0; iters <= MAXITER; iters++) {
            // COMPUTE DUALITY GAP
            double zTDDTz = F77_CALL(ddot)(&n, DTz.data(), &ione, DTz.data(), &ione);

            F77_CALL(dcopy)(&m, Dy.data(), &ione, w.data(), &ione);  // w = D*y-(mu1-mu2)
            F77_CALL(daxpy)(&m, &dminusone, mu1.data(), &ione, w.data(), &ione);
            F77_CALL(daxpy)(&m, &done, mu2.data(), &ione, w.data(), &ione);

            F77_CALL(dcopy)(&m, Dy.data(), &ione, tmp_m2.data(), &ione);  // tmp_m2 = D*y-D*D'*z
            F77_CALL(daxpy)(&m, &dminusone, DDTz.data(), &ione, tmp_m2.data(), &ione);

            F77_CALL(dcopy)(&m, w.data(), &ione, tmp_m1.data(), &ione);  // tmp_m1 = w*(D*D')^-1*w
            if (ddtf_chol) {
                F77_CALL(dpbtrs)("L", &m, &itwo, &ione, DDTF.data(), &ithree, tmp_m1.data(), &m, &info);
            } else {
                F77_CALL(dgbtrs)
                ("N", &m, &itwo, &itwo, &ione, DDTF.data(), &iseven, IPIV.data(), tmp_m1.data(), &m, &info);
            }

            double accum1 = std::accumulate(mu1.begin(), mu1.end(), 0.0);
            double accum2 = std::accumulate(mu2.begin(), mu2.end(), 0.0);
            double pobj1 = 0.5 * F77_CALL(ddot)(&m, w.data(), &ione, tmp_m1.data(), &ione) + lambda * (accum1 + accum2);
            double pobj2 = 0.5 * zTDDTz + lambda * F77_CALL(dasum)(&m, tmp_m2.data(), &ione);
            double pobj = std::min(pobj1, pobj2);
            double dobj = -0.5 * zTDDTz + F77_CALL(ddot)(&m, Dy.data(), &ione, z.data(), &ione);
            double gap = pobj - dobj;

            if (is_debug) printf("%6d %.8f %.8f %.8f\n", iters, pobj, dobj, gap);

            // STOPPING CRITERION
            if (gap <= TOL) {
                // Solved
                F77_CALL(dcopy)(&n, y, &ione, x, &ione);
                F77_CALL(daxpy)(&n, &dminusone, DTz.data(), &ione, x, &ione);
                return;
            }

            if (step >= 0.2) {
                t = std::max(2 * m * MU / gap, 1.2 * t);
            }

            // CALCULATE NEWTON STEP
            F77_CALL(dcopy)(&m, DDTz.data(), &ione, rz.data(), &ione);  // rz = D*D'*z-w
            F77_CALL(daxpy)(&m, &dminusone, w.data(), &ione, rz.data(), &ione);

            yainvx(m, +1.0 / t, f1.data(), dz.data());  // dz = r = D*y-D*D'*z+(1/t)./f1-(1/t)./f2
            yainvx(m, -1.0 / t, f2.data(), tmp_m1.data());
            F77_CALL(daxpy)(&m, &done, tmp_m1.data(), &ione, dz.data(), &ione);
            F77_CALL(daxpy)(&m, &done, tmp_m2.data(), &ione, dz.data(), &ione);

            dptr = S.data();  // S = D*D' - diag(mu1./f1-mu2./f2)
            for (int i = 0; i < m; i++) {
                *dptr++ = 6 - mu1[i] / f1[i] - mu2[i] / f2[i];
                *dptr++ = -4.0;
                *dptr++ = 1.0;
            }

            F77_CALL(dpbsv)("L", &m, &itwo, &ione, S.data(), &ithree, dz.data(), &m, &info);  // dz=S\r

            double norm2_res = F77_CALL(ddot)(&m, rz.data(), &ione, rz.data(), &ione);
            for (int i = 0; i < m; i++) {
                double tmp1 = -mu1[i] * f1[i] - (1 / t);
                double tmp2 = -mu2[i] * f2[i] - (1 / t);
                norm2_res += tmp1 * tmp1 + tmp2 * tmp2;

                dmu1[i] = -(mu1[i] + ((1 / t) + dz[i] * mu1[i]) / f1[i]);
                dmu2[i] = -(mu2[i] + ((1 / t) - dz[i] * mu2[i]) / f2[i]);
            }
            norm2_res = std::sqrt(norm2_res);

            // BACKTRACKING LINE SEARCH
            double ratio = 2;  // any number larger than 1/0.99
            for (int i = 0; i < m; i++) {
                if (dmu1[i] < 0 && -mu1[i] / dmu1[i] < ratio) ratio = -mu1[i] / dmu1[i];
                if (dmu2[i] < 0 && -mu2[i] / dmu2[i] < ratio) ratio = -mu2[i] / dmu2[i];
            }
            step = std::min(1.0, 0.99 * ratio);

            /* compute new values of z, dmu1, dmu2, f1, f2 */
            F77_CALL(daxpy)(&m, &step, dz.data(), &ione, z.data(), &ione);
            F77_CALL(daxpy)(&m, &step, dmu1.data(), &ione, mu1.data(), &ione);
            F77_CALL(daxpy)(&m, &step, dmu2.data(), &ione, mu2.data(), &ione);

            for (int lsiters = 0; lsiters < MAXLSITER; lsiters++) {  // line search loop
                int linesearch_skip = 0;
                for (int i = 0; i < m; i++) {
                    f1[i] = z[i] - lambda;
                    f2[i] = -z[i] - lambda;
                    if (f1[i] > 0 || f2[i] > 0) linesearch_skip = 1;
                }
                if (linesearch_skip != 1) {
                    DTx(m, z.data(), DTz.data());  // rz = D*D'*z-D*y-(mu1-mu2)
                    Dx(n, DTz.data(), DDTz.data());
                    F77_CALL(dcopy)(&m, DDTz.data(), &ione, rz.data(), &ione);
                    F77_CALL(daxpy)(&m, &dminusone, Dy.data(), &ione, rz.data(), &ione);
                    F77_CALL(daxpy)(&m, &done, mu1.data(), &ione, rz.data(), &ione);
                    F77_CALL(daxpy)(&m, &dminusone, mu2.data(), &ione, rz.data(), &ione);

                    // UPDATE RESIDUAL

                    // compute  norm([rz; -mu1.*f1-1/t; -mu2.*f2-1/t])
                    double norm2_newres = F77_CALL(ddot)(&m, rz.data(), &ione, rz.data(), &ione);
                    for (int i = 0; i < m; i++) {
                        double tmp1, tmp2;
                        tmp1 = -mu1[i] * f1[i] - (1 / t);
                        tmp2 = -mu2[i] * f2[i] - (1 / t);
                        norm2_newres += tmp1 * tmp1 + tmp2 * tmp2;
                    }
                    norm2_newres = sqrt(norm2_newres);

                    if (norm2_newres <= (1 - ALPHA * step) * norm2_res) break;
                }
                double diff_step = -step * (1.0 - BETA);
                F77_CALL(daxpy)(&m, &diff_step, dz.data(), &ione, z.data(), &ione);
                F77_CALL(daxpy)(&m, &diff_step, dmu1.data(), &ione, mu1.data(), &ione);
                F77_CALL(daxpy)(&m, &diff_step, dmu2.data(), &ione, mu2.data(), &ione);
                step *= BETA;
            }
        }

        // Maxiter exceeded
        F77_CALL(dcopy)(&n, y, &ione, x, &ione);
        F77_CALL(daxpy)(&n, &dminusone, DTz.data(), &ione, x, &ione);
    }

    void calc() {
        lambda_max = l1tf_lambdamax(m_y.size(), m_y.data());
        // printf("calc lambda (lambda_max) = %e (%e)\n\n", lambda, lambda_max);
        l1tf(m_y.size(), m_y.data(), m_x.data());
    }

    void set_data(const std::vector<double> &y_) {
        m_y = y_;
        m_x.resize(m_y.size(), 0);
    }

    bool is_debug{false};
    std::vector<double> m_y, m_x;
    double lambda{0.1}, lambda_max{0};
};

}  // namespace ornate

#endif
