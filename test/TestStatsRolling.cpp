#include <math_dummy.h>
#include <math_stats_rolling_rb_range.h>
#include <math_utils.h>
#include <iostream>
#include "catch.hpp"
#include "math_stats.h"

using namespace std;
using namespace ornate;

void test_variance_rolling(const vector<double>& _data, int window) {
    rolling_variance_rb_range rvrr(1);
    rolling_variance_rb_range nrvrr(1);

    for (int round = 0; round < 2; ++round) {
        rvrr.init();
        nrvrr.init();

        ornate::variance_rolling<> vr(window);
        rolling_variance_rb rvr(window);
        rolling_data_container<> container(window, 1);
        vector<double> row(1, 0);

        for (int i = 0; i < 6; ++i) {
            int to = i + 1;
            int from = to - window;
            if (from < 0) from = 0;
            double std_naive = ornate::std(_data.data() + from, to - from);
            double s = sqrtf(vr(_data[i]));
            rvr(_data[i]);
            row[0] = _data[i];
            container.push(row);
            rvrr(container.get_old_row(), container.get_new_row(), row.data());

            REQUIRE(FloatEqual(sqrtf(vr.variance), std_naive));
            REQUIRE(FloatEqual(sqrtf(rvr.variance), std_naive));
            REQUIRE(FloatEqual(sqrtf(row[0]), std_naive));

            int real_window = window;
            if (container.m_count < window) {
                real_window = container.m_count;
            }
            nrvrr.init();
            for (int ts_idx = 0; ts_idx < real_window; ++ts_idx) {
                nrvrr.full_single(ts_idx, container.get_row_by_idx(ts_idx), row.data());
            }
            nrvrr.final_result(row.data());
            REQUIRE(FloatEqual(sqrtf(row[0]), std_naive));
        }
    }
}

void test_covariance_rolling(const vector<double>& _data1, const vector<double>& _data2, int window) {
    rolling_cov_rb_range rcrr(1);
    rolling_cov_rb_range nrcrr(1);

    for (int round = 0; round < 2; ++round) {
        rcrr.init();
        nrcrr.init();

        ornate::covariance_rolling<> cvr(window);
        rolling_cov_rb rcr(window);
        rolling_data_container<> container(window, 1);
        rolling_data_container<> container2(window, 1);
        vector<double> row(1, 0);
        vector<double> row2(1, 0);

        for (int i = 0; i < 6; ++i) {
            int to = i + 1;
            int from = to - window;
            if (from < 0) from = 0;
            double cov_naive = ornate::cov(_data1.data() + from, _data2.data() + from, to - from);
            double s = cvr(_data1[i], _data2[i]);
            rcr(_data1[i], _data2[i]);
            row[0] = _data1[i];
            container.push(row);
            row2[0] = _data2[i];
            container2.push(row2);
            rcrr(container.get_old_row(), container2.get_old_row(), container.get_new_row(), container2.get_new_row(),
                 row.data());

            REQUIRE(FloatEqual(cvr.covariance, cov_naive));
            REQUIRE(FloatEqual(rcr.covariance, cov_naive));
            REQUIRE(FloatEqual(row[0], cov_naive));

            int real_window = window;
            if (container.m_count < window) {
                real_window = container.m_count;
            }

            nrcrr.init();
            for (int ts_idx = 0; ts_idx < real_window; ++ts_idx) {
                nrcrr.full_single(ts_idx, container.get_row_by_idx(ts_idx), container2.get_row_by_idx(ts_idx),
                                  row.data());
            }
            nrcrr.final_result(row.data());
            REQUIRE(FloatEqual(row[0], cov_naive));
        }
    }
}

void test_corr_rolling(const vector<double>& _data1, const vector<double>& _data2, int window) {
    rolling_data_container<> container(window, 1);
    rolling_data_container<> container2(window, 1);
    vector<double> row(1, 0);
    vector<double> row2(1, 0);
    rolling_corr_rb_range rcrr(1);
    rolling_corr_rb_range nrcrr(1);

    for (int round = 0; round < 2; ++round) {
        container.clear();
        container2.clear();
        rcrr.init();
        nrcrr.init();

        ornate::corr_rolling cvr(window);
        rolling_corr_rb rcr(window);
        for (int i = 0; i < 6; ++i) {
            int to = i + 1;
            int from = to - window;
            if (from < 0) from = 0;
            double corr_naive = ornate::corr(_data1.data() + from, _data2.data() + from, to - from);
            double s = cvr(_data1[i], _data2[i]);
            rcr(_data1[i], _data2[i]);

            row[0] = _data1[i];
            container.push(row);
            row2[0] = _data2[i];
            container2.push(row2);
            rcrr(container.get_old_row(), container2.get_old_row(), container.get_new_row(), container2.get_new_row(),
                 row.data());

            REQUIRE(FloatEqual(s, corr_naive));
            REQUIRE(FloatEqual(rcr.corr, corr_naive));
            REQUIRE(FloatEqual(row[0], corr_naive));

            int real_window = window;
            if (container.m_count < window) {
                real_window = container.m_count;
            }

            nrcrr.init();
            for (int ts_idx = 0; ts_idx < real_window; ++ts_idx) {
                nrcrr.full_single(ts_idx, container.get_row_by_idx(ts_idx), container2.get_row_by_idx(ts_idx),
                                  row.data());
            }
            nrcrr.final_result(row.data());
            REQUIRE(FloatEqual(row[0], corr_naive));
        }
    }
}

void test_skew_by_window(const vector<double>& x_, int window) {
    rolling_data_container<> container(window, 1);
    vector<double> row(1, 0);
    rolling_skew_rb_range rsrr(1);
    rolling_skew_rb_range nrsrr(1);

    for (int round = 0; round < 2; ++round) {
        container.clear();
        rsrr.init();
        nrsrr.init();

        rolling_skew_rb sr(window);
        vector<double> y;

        double ret = 0;
        for (double d : x_) {
            add_window_vector(y, window, d);
            ret = sr(d);

            row[0] = d;
            container.push(row);
            rsrr(container.get_old_row(), container.get_new_row(), row.data());

            double expected = math_skew(y);
            REQUIRE(FloatEqual(ret, expected));
            REQUIRE(FloatEqual(row[0], expected));

            int real_window = window;
            if (container.m_count < window) {
                real_window = container.m_count;
            }

            nrsrr.init();
            for (int ts_idx = 0; ts_idx < real_window; ++ts_idx) {
                nrsrr.full_single(ts_idx, container.get_row_by_idx(ts_idx), row.data());
            }
            nrsrr.final_result(row.data());
            REQUIRE(FloatEqual(row[0], expected));
        }
    }
}

void test_kurtosis_by_window(const vector<double>& x_, int window) {
    rolling_data_container<> container(window, 1);
    vector<double> row(1, 0);
    rolling_kurtosis_rb_range rkrr(1);
    rolling_kurtosis_rb_range nrkrr(1);

    for (int round = 0; round < 2; ++round) {
        container.clear();
        rkrr.init();
        nrkrr.init();

        rolling_kurtosis_rb sr(window);
        vector<double> y;

        double ret = 0;
        for (double d : x_) {
            add_window_vector(y, window, d);
            ret = sr(d);

            row[0] = d;
            container.push(row);
            rkrr(container.get_old_row(), container.get_new_row(), row.data());

            double expected = math_kurtosis(y);

            REQUIRE(FloatEqual(ret, expected));
            REQUIRE(FloatEqual(row[0], expected));

            int real_window = window;
            if (container.m_count < window) {
                real_window = container.m_count;
            }

            nrkrr.init();
            for (int ts_idx = 0; ts_idx < real_window; ++ts_idx) {
                nrkrr.full_single(ts_idx, container.get_row_by_idx(ts_idx), row.data());
            }
            nrkrr.final_result(row.data());
            REQUIRE(FloatEqual(ret, row[0]));
        }
    }
}

void test_decay_by_window(const vector<double>& x_, int window) {
    rolling_data_container<> container(window, 1);
    vector<double> row(1, 0);
    rolling_decay_rb_range rdrr(1);
    rdrr.set_row_size(window);
    rolling_decay_rb_range nrdrr(1);
    nrdrr.set_row_size(window);

    for (int round = 0; round < 2; ++round) {
        container.clear();
        rdrr.init();
        nrdrr.init();

        vector<double> y;
        for (double d : x_) {
            add_window_vector(y, window, d);

            row[0] = d;
            container.push(row);
            rdrr(container.get_old_row(), container.get_new_row(), row.data());

            double expected = dummy_decay(y);
            REQUIRE(FloatEqual(row[0], expected));

            int real_window = window;
            if (container.m_count < window) {
                real_window = container.m_count;
            }

            nrdrr.init();
            for (int ts_idx = 0; ts_idx < real_window; ++ts_idx) {
                nrdrr.full_single(ts_idx, real_window, container.get_row_by_idx(ts_idx), row.data());
            }
            nrdrr.final_result(row.data());
            REQUIRE(FloatEqual(expected, row[0]));
        }
    }
}

void test_rank_by_window(const vector<double>& x_, int window) {
    rolling_data_container<> container(window, 2);
    vector<double> row(2, 0);
    rolling_rank_count_rb_range<double> rrrr(2);
    rrrr.set_row_size(window);

    for (int round = 0; round < 2; ++round) {
        container.clear();
        rrrr.init();

        rolling_rank_count_rb<double> rrrb(window);
        vector<double> y;
        double ret = 0;
        for (double d : x_) {
            add_window_vector(y, window, d);
            ret = rrrb(d);

            row[0] = d;
            row[1] = d;
            container.push(row);
            rrrr(container.get_old_row(), container.get_new_row(), row.data());

            double expected = ornate::rank_last(y.data(), y.size() - 1, window);

            REQUIRE(FloatEqual(ret, expected));
            REQUIRE(FloatEqual(row[0], expected));
            REQUIRE(FloatEqual(row[1], expected));
        }
    }
}

void test_regression_by_window(const vector<double>& x_, const vector<double>& y_, int window) {
    rolling_data_container<> container(window, 2);
    rolling_data_container<> container2(window, 2);
    vector<double> row(2, 0);
    vector<double> row2(2, 0);
    rolling_regression2_rb_range rrrr(2);
    rrrr.set_row_size(window);
    rolling_regression2_rb_range nrrrr(2);
    nrrrr.set_row_size(window);

    for (int round = 0; round < 2; ++round) {
        container.clear();
        container2.clear();
        rrrr.init();
        nrrrr.init();

        regression2_rolling rr(window);
        vector<double> _x, _y;
        rolling_regression2_rb rrrb(window);

        double ret = 0;
        for (size_t i = 0; i < x_.size(); ++i) {
            add_window_vector(_x, window, x_[i]);
            add_window_vector(_y, window, y_[i]);

            double a, b, R;
            ornate::regression(_y, _x, &a, &b, &R);

            rr(y_[i], x_[i]);
            rrrb(y_[i], x_[i]);

            REQUIRE(FloatEqual(rr.a, a));
            REQUIRE(FloatEqual(rr.b, b));
            REQUIRE(FloatEqual(rrrb.a, a));
            REQUIRE(FloatEqual(rrrb.b, b));

            row[0] = y_[i];
            row[1] = y_[i];
            container.push(row);
            row2[0] = x_[i];
            row2[1] = x_[i];
            container2.push(row2);
            rrrr(container.get_old_row(), container2.get_old_row(), container.get_new_row(), container2.get_new_row());
            rrrr.get_coefficients(row.data(), row2.data());
            REQUIRE(FloatEqual(row[0], a));
            REQUIRE(FloatEqual(row[1], a));
            REQUIRE(FloatEqual(row2[0], b));
            REQUIRE(FloatEqual(row2[1], b));
            rrrr.get_fitted(container2.get_new_row(), row2.data());
            REQUIRE(FloatEqual(row2[0], a + b * x_[i]));
            rrrr.get_residual(container.get_new_row(), container2.get_new_row(), row2.data());
            REQUIRE(FloatEqual(row2[0], y_[i] - (a + b * x_[i])));

            int real_window = window;
            if (container.m_count < window) {
                real_window = container.m_count;
            }

            rrrr.r2_pre_calc();
            for (int ts_idx = 0; ts_idx < real_window; ++ts_idx) {
                rrrr.r2_single(ts_idx, container.get_row_by_idx(ts_idx), container2.get_row_by_idx(ts_idx));
            }
            rrrr.get_r2(row2.data());
            double expected_r2 = dummy_r2(_x, _y, a, b, real_window);
            REQUIRE(FloatEqual(row2[0], expected_r2));

            {
                nrrrr.init();

                for (int ts_idx = 0; ts_idx < real_window; ++ts_idx) {
                    nrrrr.full_single(ts_idx, container.get_row_by_idx(ts_idx), container2.get_row_by_idx(ts_idx));
                }
                nrrrr.get_coefficients(row.data(), row2.data());
                REQUIRE(FloatEqual(a, row[0]));
                REQUIRE(FloatEqual(b, row2[0]));

                nrrrr.get_fitted(container2.get_new_row(), row2.data());
                REQUIRE(FloatEqual(row2[0], a + b * x_[i]));

                nrrrr.get_residual(container.get_new_row(), container2.get_new_row(), row2.data());
                REQUIRE(FloatEqual(row2[0], y_[i] - (a + b * x_[i])));

                nrrrr.r2_pre_calc();
                for (int ts_idx = 0; ts_idx < real_window; ++ts_idx) {
                    nrrrr.r2_single(ts_idx, container.get_row_by_idx(ts_idx), container2.get_row_by_idx(ts_idx));
                }
                nrrrr.get_r2(row2.data());
                expected_r2 = dummy_r2(_x, _y, a, b, real_window);
                REQUIRE(FloatEqual(row2[0], expected_r2));
            }
        }
    }
}

void test_ema_hl_by_window(const vector<double>& x_, int window) {
    rolling_data_container<> container(window, 2);
    vector<double> row(2, 0);
    rolling_ema_hl_rb_range rqrr(2);
    rqrr.set_row_size(window);
    rqrr.set_param("hl", "4");

    rolling_ema_hl_rb_range nrerr(2);
    nrerr.set_row_size(window);
    nrerr.set_param("hl", "4");

    for (int round = 0; round < 2; ++round) {
        container.clear();
        rqrr.init();
        nrerr.init();

        rolling_ema_hl_rb rerb(window, 4);
        vector<double> y;

        double ret = 0;
        for (double d : x_) {
            add_window_vector(y, window, d);
            ret = rerb(d);

            row[0] = d;
            row[1] = d;
            container.push(row);
            rqrr(container.get_old_row(), container.get_new_row(), row.data());

            vector<double> tmp = y;
            double expected = ornate::ema_hl(tmp, window, window - 1, 4);

            REQUIRE(FloatEqual(ret, expected));
            REQUIRE(FloatEqual(row[0], expected));
            REQUIRE(FloatEqual(row[1], expected));

            int real_window = window;
            if (container.m_count < window) {
                real_window = container.m_count;
            }

            nrerr.init();
            for (int ts_idx = 0; ts_idx < real_window; ++ts_idx) {
                nrerr.full_single(ts_idx, container.get_row_by_idx(ts_idx), row.data());
            }
            nrerr.final_result(row.data());
            REQUIRE(FloatEqual(expected, row[0]));
        }
    }
}

void test_regression3_by_window(const vector<double>& x1_, const vector<double>& x2_, const vector<double>& y_,
                                int window) {
    rolling_data_container<> container(window, 2);
    rolling_data_container<> container2(window, 2);
    rolling_data_container<> container3(window, 2);
    vector<double> row(2, 0);
    vector<double> row2(2, 0);
    vector<double> row3(2, 0);
    rolling_regression3_rb_range rrrr(2);
    rrrr.set_row_size(window);
    rolling_regression3_rb_range nrrrr(2);
    nrrrr.set_row_size(window);

    for (int round = 0; round < 2; ++round) {
        container.clear();
        container2.clear();
        container3.clear();
        rrrr.init();
        nrrrr.init();

        vector<double> _x1, _x2, _y;
        rolling_regression3_rb rrrb(window);

        double ret = 0;
        for (size_t i = 0; i < x1_.size(); ++i) {
            add_window_vector(_x1, window, x1_[i]);
            add_window_vector(_x2, window, x2_[i]);
            add_window_vector(_y, window, y_[i]);

            double b0, b1, b2;
            ornate::regression3(_y, _x1, _x2, &b0, &b1, &b2);

            rrrb(y_[i], x1_[i], x2_[i]);
            REQUIRE(FloatEqual(rrrb.b0, b0));
            REQUIRE(FloatEqual(rrrb.b1, b1));
            REQUIRE(FloatEqual(rrrb.b2, b2));

            row[0] = y_[i];
            row[1] = y_[i];
            container.push(row);
            row2[0] = x1_[i];
            row2[1] = x1_[i];
            container2.push(row2);
            row3[0] = x2_[i];
            row3[1] = x2_[i];
            container3.push(row3);
            rrrr(container.get_old_row(), container2.get_old_row(), container3.get_old_row(), container.get_new_row(),
                 container2.get_new_row(), container3.get_new_row());
            rrrr.get_coefficients(row.data(), row2.data(), row3.data());
            REQUIRE(FloatEqual(row[0], b0));
            REQUIRE(FloatEqual(row[1], b0));
            REQUIRE(FloatEqual(row2[0], b1));
            REQUIRE(FloatEqual(row2[1], b1));
            REQUIRE(FloatEqual(row3[0], b2));
            REQUIRE(FloatEqual(row3[1], b2));

            rrrr.get_fitted(container2.get_new_row(), container3.get_new_row(), row2.data());
            REQUIRE(FloatEqual(row2[0], b0 + b1 * x1_[i] + b2 * x2_[i]));
            rrrr.get_residual(container.get_new_row(), container2.get_new_row(), container3.get_new_row(), row2.data());
            REQUIRE(FloatEqual(row2[0], y_[i] - (b0 + b1 * x1_[i] + b2 * x2_[i])));

            int real_window = window;
            if (container.m_count < window) {
                real_window = container.m_count;
            }

            rrrr.r2_pre_calc();
            for (int ts_idx = 0; ts_idx < real_window; ++ts_idx) {
                rrrr.r2_single(ts_idx, container.get_row_by_idx(ts_idx), container2.get_row_by_idx(ts_idx),
                               container3.get_row_by_idx(ts_idx));
            }
            rrrr.get_r2(row2.data());
            double expected_r2 = dummy_r2(_x1, _x2, _y, b0, b1, b2, real_window);
            REQUIRE(FloatEqual(row2[0], expected_r2));

            {
                nrrrr.init();

                for (int ts_idx = 0; ts_idx < real_window; ++ts_idx) {
                    nrrrr.full_single(ts_idx, container.get_row_by_idx(ts_idx), container2.get_row_by_idx(ts_idx),
                                      container3.get_row_by_idx(ts_idx));
                }
                nrrrr.get_coefficients(row.data(), row2.data(), row3.data());
                REQUIRE(FloatEqual(b0, row[0]));
                REQUIRE(FloatEqual(b1, row2[0]));
                REQUIRE(FloatEqual(b2, row3[0]));

                nrrrr.get_fitted(container2.get_new_row(), container3.get_new_row(), row2.data());
                REQUIRE(FloatEqual(row2[0], b0 + b1 * x1_[i] + b2 * x2_[i]));
                nrrrr.get_residual(container.get_new_row(), container2.get_new_row(), container3.get_new_row(),
                                   row2.data());
                REQUIRE(FloatEqual(row2[0], y_[i] - (b0 + b1 * x1_[i] + b2 * x2_[i])));

                nrrrr.r2_pre_calc();
                for (int ts_idx = 0; ts_idx < real_window; ++ts_idx) {
                    nrrrr.r2_single(ts_idx, container.get_row_by_idx(ts_idx), container2.get_row_by_idx(ts_idx),
                                    container3.get_row_by_idx(ts_idx));
                }
                nrrrr.get_r2(row2.data());
                expected_r2 = dummy_r2(_x1, _x2, _y, b0, b1, b2, real_window);
                REQUIRE(FloatEqual(row2[0], expected_r2));
            }
        }
    }
}

void test_quantile_by_window(const vector<double>& x_, int window, double percent) {
    rolling_data_container<> container(window, 2);
    vector<double> row(2, 0);
    rolling_quantile_rb_range<double> rqrr(2, percent);
    rqrr.set_row_size(window);

    for (int round = 0; round < 2; ++round) {
        container.clear();
        rqrr.init();

        rolling_quantile_rb<double> rrrb(window, percent);
        rolling_weighted_quantile_rb<double> rwqrb(window, percent);
        vector<double> y;
        double ret = 0, ret1 = 0;
        for (double d : x_) {
            add_window_vector(y, window, d);
            ret = rrrb(d);
            ret1 = rwqrb(d, 1.0);

            row[0] = d;
            row[1] = d;
            container.push(row);
            rqrr(container.get_old_row(), container.get_new_row(), row.data());

            vector<double> tmp = y;
            double expected = ornate::quantile(tmp, percent);

            REQUIRE(FloatEqual(ret, expected));
            REQUIRE(FloatEqual(ret1, expected));
            REQUIRE(FloatEqual(row[0], expected));
            REQUIRE(FloatEqual(row[1], expected));
        }
    }
}

void test_ema_hl_pp_by_window(const vector<double>& x_, int window) {
    rolling_data_container<> container(window, 1);
    vector<double> row(1, 0);
    rolling_ema_hl_rb_range rqrr(1);
    rqrr.set_row_size(window);
    rqrr.set_param("param1", "4");
    rqrr.set_param("param2", "8");

    rolling_ema_hl_rb_range nrerr(1);
    nrerr.set_row_size(window);
    nrerr.set_param("param1", "4");
    nrerr.set_param("param2", "8");

    for (int round = 0; round < 2; ++round) {
        container.clear();
        rqrr.init();
        nrerr.init();

        double ret = 0;
        for (double d : x_) {
            row[0] = d;
            container.push(row);
            rqrr(container.get_old_row(), container.get_new_row(), row.data());

            double result = row[0];

            int real_window = window;
            if (container.m_count < window) {
                real_window = container.m_count;
            }

            nrerr.init();
            for (int ts_idx = 0; ts_idx < real_window; ++ts_idx) {
                nrerr.full_single(ts_idx, container.get_row_by_idx(ts_idx), row.data());
            }
            nrerr.final_result(row.data());
            REQUIRE(FloatEqual(result, row[0]));
        }
    }
}

void test_ols_by_window(const vector<double>& x_, const vector<double>& y_, int window) {
    rolling_data_container<> container(window, 2);
    rolling_data_container<> container2(window, 2);
    vector<double> row(2, 0);
    vector<double> row2(2, 0);
    rolling_ols2_rb_range rrrr(2);
    rrrr.set_row_size(window);
    rolling_ols2_rb_range nrrrr(2);
    nrrrr.set_row_size(window);

    for (int round = 0; round < 2; ++round) {
        container.clear();
        container2.clear();
        rrrr.init();
        nrrrr.init();

        ols2_rolling rr(window);
        vector<double> _x, _y;
        rolling_ols2_rb rrrb(window);
        double ret = 0;
        for (size_t i = 0; i < x_.size(); ++i) {
            add_window_vector(_x, window, x_[i]);
            add_window_vector(_y, window, y_[i]);

            double b = ornate::ols(_y, _x);

            rr(y_[i], x_[i]);
            rrrb(y_[i], x_[i]);

            REQUIRE(FloatEqual(rr.b, b));
            REQUIRE(FloatEqual(rrrb.b, b));

            row[0] = y_[i];
            row[1] = y_[i];
            container.push(row);
            row2[0] = x_[i];
            row2[1] = x_[i];
            container2.push(row2);
            rrrr(container.get_old_row(), container2.get_old_row(), container.get_new_row(), container2.get_new_row());
            rrrr.get_coefficients(row2.data());
            REQUIRE(FloatEqual(row2[0], b));
            REQUIRE(FloatEqual(row2[1], b));
            rrrr.get_fitted(container2.get_new_row(), row2.data());
            REQUIRE(FloatEqual(row2[0], b * x_[i]));
            rrrr.get_residual(container.get_new_row(), container2.get_new_row(), row2.data());
            REQUIRE(FloatEqual(row2[0], y_[i] - (b * x_[i])));

            int real_window = window;
            if (container.m_count < window) {
                real_window = container.m_count;
            }

            rrrr.r2_pre_calc();
            for (int ts_idx = 0; ts_idx < real_window; ++ts_idx) {
                rrrr.r2_single(ts_idx, container.get_row_by_idx(ts_idx), container2.get_row_by_idx(ts_idx));
            }
            rrrr.get_r2(row2.data());
            double expected_r2 = dummy_r2_no_slope(_x, _y, b, real_window);
            REQUIRE(FloatEqual(row2[0], expected_r2));

            {
                nrrrr.init();

                for (int ts_idx = 0; ts_idx < real_window; ++ts_idx) {
                    nrrrr.full_single(ts_idx, container.get_row_by_idx(ts_idx), container2.get_row_by_idx(ts_idx));
                }
                nrrrr.get_coefficients(row2.data());
                REQUIRE(FloatEqual(b, row2[0]));

                nrrrr.get_fitted(container2.get_new_row(), row2.data());
                REQUIRE(FloatEqual(row2[0], b * x_[i]));

                nrrrr.get_residual(container.get_new_row(), container2.get_new_row(), row2.data());
                REQUIRE(FloatEqual(row2[0], y_[i] - (b * x_[i])));

                nrrrr.r2_pre_calc();
                for (int ts_idx = 0; ts_idx < real_window; ++ts_idx) {
                    nrrrr.r2_single(ts_idx, container.get_row_by_idx(ts_idx), container2.get_row_by_idx(ts_idx));
                }
                nrrrr.get_r2(row2.data());
                expected_r2 = dummy_r2_no_slope(_x, _y, b, real_window);
                REQUIRE(FloatEqual(row2[0], expected_r2));
            }
        }
    }
}

void test_ols3_by_window(const vector<double>& x1_, const vector<double>& x2_, const vector<double>& y_, int window) {
    rolling_data_container<> container(window, 2);
    rolling_data_container<> container2(window, 2);
    rolling_data_container<> container3(window, 2);
    vector<double> row(2, 0);
    vector<double> row2(2, 0);
    vector<double> row3(2, 0);
    rolling_ols3_rb_range rrrr(2);
    rrrr.set_row_size(window);
    rolling_ols3_rb_range nrrrr(2);
    nrrrr.set_row_size(window);

    for (int round = 0; round < 2; ++round) {
        container.clear();
        container2.clear();
        container3.clear();
        rrrr.init();
        nrrrr.init();

        vector<double> _x1, _x2, _y;
        rolling_ols3_rb rrrb(window);

        double ret = 0;
        for (size_t i = 0; i < x1_.size(); ++i) {
            add_window_vector(_x1, window, x1_[i]);
            add_window_vector(_x2, window, x2_[i]);
            add_window_vector(_y, window, y_[i]);

            double b1, b2;
            ornate::ols(_y, _x1, _x2, &b1, &b2);

            rrrb(y_[i], x1_[i], x2_[i]);

            REQUIRE(FloatEqual(rrrb.b1, b1));
            REQUIRE(FloatEqual(rrrb.b2, b2));

            row[0] = y_[i];
            row[1] = y_[i];
            container.push(row);
            row2[0] = x1_[i];
            row2[1] = x1_[i];
            container2.push(row2);
            row3[0] = x2_[i];
            row3[1] = x2_[i];
            container3.push(row3);
            rrrr(container.get_old_row(), container2.get_old_row(), container3.get_old_row(), container.get_new_row(),
                 container2.get_new_row(), container3.get_new_row());
            rrrr.get_coefficients(row2.data(), row3.data());
            REQUIRE(FloatEqual(row2[0], b1));
            REQUIRE(FloatEqual(row2[1], b1));
            REQUIRE(FloatEqual(row3[0], b2));
            REQUIRE(FloatEqual(row3[1], b2));

            rrrr.get_fitted(container2.get_new_row(), container3.get_new_row(), row2.data());
            REQUIRE(FloatEqual(row2[0], b1 * x1_[i] + b2 * x2_[i]));
            rrrr.get_residual(container.get_new_row(), container2.get_new_row(), container3.get_new_row(), row2.data());
            REQUIRE(FloatEqual(row2[0], y_[i] - (b1 * x1_[i] + b2 * x2_[i])));

            int real_window = window;
            if (container.m_count < window) {
                real_window = container.m_count;
            }

            {
                rrrr.r2_pre_calc();
                for (int ts_idx = 0; ts_idx < real_window; ++ts_idx) {
                    rrrr.r2_single(ts_idx, container.get_row_by_idx(ts_idx), container2.get_row_by_idx(ts_idx),
                                   container3.get_row_by_idx(ts_idx));
                }
                rrrr.get_r2(row2.data());
                double expected_r2 = dummy_r2_no_slope(_x1, _x2, _y, b1, b2, real_window);
                REQUIRE(FloatEqual(row2[0], expected_r2));
            }

            {
                nrrrr.init();

                for (int ts_idx = 0; ts_idx < real_window; ++ts_idx) {
                    nrrrr.full_single(ts_idx, container.get_row_by_idx(ts_idx), container2.get_row_by_idx(ts_idx),
                                      container3.get_row_by_idx(ts_idx));
                }

                nrrrr.get_coefficients(row2.data(), row3.data());
                REQUIRE(FloatEqual(b1, row2[0]));
                REQUIRE(FloatEqual(b2, row3[0]));

                nrrrr.get_fitted(container2.get_new_row(), container3.get_new_row(), row2.data());
                REQUIRE(FloatEqual(row2[0], b1 * x1_[i] + b2 * x2_[i]));
                nrrrr.get_residual(container.get_new_row(), container2.get_new_row(), container3.get_new_row(),
                                   row2.data());
                REQUIRE(FloatEqual(row2[0], y_[i] - (b1 * x1_[i] + b2 * x2_[i])));

                nrrrr.r2_pre_calc();
                for (int ts_idx = 0; ts_idx < real_window; ++ts_idx) {
                    nrrrr.r2_single(ts_idx, container.get_row_by_idx(ts_idx), container2.get_row_by_idx(ts_idx),
                                    container3.get_row_by_idx(ts_idx));
                }
                nrrrr.get_r2(row2.data());
                double expected_r2 = dummy_r2_no_slope(_x1, _x2, _y, b1, b2, real_window);
                REQUIRE(FloatEqual(row2[0], expected_r2));
            }
        }
    }
}

void test_slope_no_intercept_by_window(const vector<double>& y_, int window) {
    rolling_data_container<> container(window, 2);
    vector<double> row(2, 0);
    slope_no_intercept_rb_range nrrrr(2);
    nrrrr.set_row_size(window);

    for (int round = 0; round < 2; ++round) {
        container.clear();
        nrrrr.init();

        slope_no_intercept_rolling rr(window);
        vector<double> _y;

        double ret = 0;
        for (double i : y_) {
            add_window_vector(_y, window, i);

            double b = ornate::slope_no_intercept(_y);
            double b1 = rr(i);
            REQUIRE(FloatEqual(b1, b));

            row[0] = i;
            row[1] = i;
            container.push(row);

            int real_window = window;
            if (container.m_count < window) {
                real_window = container.m_count;
            }

            {
                nrrrr.init();
                for (int ts_idx = 0; ts_idx < real_window; ++ts_idx) {
                    nrrrr.full_single(ts_idx, real_window, container.get_row_by_idx(ts_idx), row.data());
                }
                nrrrr.final_result(row.data());
                REQUIRE(FloatEqual(b, row[0]));
            }
        }
    }
}

void test_slope_by_window(const vector<double>& y_, int window) {
    rolling_data_container<> container(window, 2);
    vector<double> row(2, 0);
    slope_rb_range nrrrr(2);
    nrrrr.set_row_size(window);

    for (int round = 0; round < 2; ++round) {
        container.clear();
        nrrrr.init();

        slope_rolling rr(window);
        vector<double> _y;
        double ret = 0;
        for (double i : y_) {
            add_window_vector(_y, window, i);

            double b = ornate::ts_slope(_y);

            double b1 = rr(i);
            REQUIRE(FloatEqual(b1, b));

            row[0] = i;
            row[1] = i;
            container.push(row);

            int real_window = window;
            if (container.m_count < window) {
                real_window = container.m_count;
            }

            {
                nrrrr.init();
                for (int ts_idx = 0; ts_idx < real_window; ++ts_idx) {
                    nrrrr.full_single(ts_idx, real_window, container.get_row_by_idx(ts_idx), row.data());
                }
                nrrrr.final_result(row.data());
                REQUIRE(FloatEqual(b, row[0]));
            }
        }
    }
}

void test_sharpe_by_window(const vector<double>& y_, int window) {
    rolling_data_container<> container(window, 2);
    vector<double> row(2, 0);
    rolling_sharpe_rb_range rsrr(2);
    rsrr.set_row_size(window);
    rolling_sharpe_rb_range nrsrr(2);
    nrsrr.set_row_size(window);

    for (int round = 0; round < 2; ++round) {
        container.clear();
        rsrr.init();
        nrsrr.init();

        vector<double> _y;
        double ret = 0;
        for (double i : y_) {
            add_window_vector(_y, window, i);

            double b = ornate::ts_sharpe(_y);

            row[0] = i;
            row[1] = i;
            container.push(row);
            rsrr(container.get_old_row(), container.get_new_row(), row.data());
            REQUIRE(FloatEqual(row[0], b));
            REQUIRE(FloatEqual(row[1], b));

            int real_window = window;
            if (container.m_count < window) {
                real_window = container.m_count;
            }

            {
                nrsrr.init();
                for (int ts_idx = 0; ts_idx < real_window; ++ts_idx) {
                    nrsrr.full_single(ts_idx, container.get_row_by_idx(ts_idx), row.data());
                }
                nrsrr.final_result(row.data());
                REQUIRE(FloatEqual(b, row[0]));
            }
        }
    }
}

void test_scale_by_window(const vector<double>& y_, int window) {
    rolling_data_container<> container(window, 2);
    vector<double> row(2, 0);
    rolling_scale_rb_range rsrr(2);
    rsrr.set_row_size(window);
    rolling_scale_rb_range nrsrr(2);
    nrsrr.set_row_size(window);

    for (int round = 0; round < 2; ++round) {
        container.clear();
        rsrr.init();
        nrsrr.init();

        vector<double> _y;
        double ret = 0;
        for (double i : y_) {
            add_window_vector(_y, window, i);

            double b = ornate::ts_scale(_y);

            row[0] = i;
            row[1] = i;
            container.push(row);
            rsrr(container.get_old_row(), container.get_new_row(), row.data());
            REQUIRE(FloatEqual(row[0], b));
            REQUIRE(FloatEqual(row[1], b));

            int real_window = window;
            if (container.m_count < window) {
                real_window = container.m_count;
            }

            {
                nrsrr.init();
                for (int ts_idx = 0; ts_idx < real_window; ++ts_idx) {
                    nrsrr.full_single(ts_idx, container.get_row_by_idx(ts_idx), row.data());
                }
                nrsrr.final_result(row.data());
                REQUIRE(FloatEqual(b, row[0]));
            }
        }
    }
}

void test_rank2_by_window(const vector<double>& x_, const vector<double>& x1_, int window) {
    rolling_data_container<> container(window, 2);
    vector<double> row(2, 0);
    rolling_data_container<> container1(window, 2);
    vector<double> row1(2, 0);
    rolling_rank2_count_rb_range<double> rrrr(2);
    rrrr.set_row_size(window);

    for (int round = 0; round < 2; ++round) {
        container.clear();
        container1.clear();
        rrrr.init();

        vector<double> y;
        double ret = 0;
        for (size_t i = 0; i < x_.size(); ++i) {
            add_window_vector(y, window, x1_[i]);

            row[0] = x_[i];
            row[1] = x_[i];
            container.push(row);
            row1[0] = x1_[i];
            row1[1] = x1_[i];
            container1.push(row1);
            rrrr(container.get_old_row(), container1.get_old_row(), container.get_new_row(), container1.get_new_row(),
                 row.data());

            double expected = ornate::rank2_last(y.data(), x_[i], y.size() - 1, window);

            REQUIRE(FloatEqual(row[0], expected));
            REQUIRE(FloatEqual(row[1], expected));
        }
    }
}

void test_dcor_rolling(const vector<double>& _data1, const vector<double>& _data2, int window) {
    rolling_data_container<> container(window, 1);
    rolling_data_container<> container2(window, 1);
    vector<double> row(1, 0);
    vector<double> row2(1, 0);
    rolling_dcor_rb_range nrdrr(1);
    nrdrr.set_row_size(window);
    rolling_dcor_rb_range rdrr(1);
    rdrr.set_row_size(window);

    for (int round = 0; round < 2; ++round) {
        container.clear();
        container2.clear();
        nrdrr.init();
        rdrr.init();

        vector<double> x_, y;
        x_.push_back(_data1[0]);
        y.push_back(_data2[0]);
        for (size_t i = 1; i < _data1.size(); ++i) {
            double d = _data1[i] - _data1[i - 1];
            add_window_vector(x_, window + 1, _data1[i]);
            double d1 = _data2[i] - _data2[i - 1];
            add_window_vector(y, window + 1, _data2[i]);

            row[0] = d;
            container.push(row);
            row2[0] = d1;
            container2.push(row2);
            rdrr(container.get_old_row(), container2.get_old_row(), container.get_new_row(), container2.get_new_row(),
                 row.data());

            double corr_naive = dummy_dcor(x_, y);
            REQUIRE(FloatEqual(row[0], corr_naive));

            int real_window = window;
            if (container.m_count < window) {
                real_window = container.m_count;
            }

            {
                nrdrr.init();
                for (int ts_idx = 0; ts_idx < real_window; ++ts_idx) {
                    nrdrr.full_single(ts_idx, container.get_row_by_idx(ts_idx), container2.get_row_by_idx(ts_idx),
                                      row.data());
                }
                nrdrr.final_result(row.data());
                REQUIRE(FloatEqual(row[0], corr_naive));
            }
        }
    }
}

void test_tscross_rolling(const vector<double>& _data1, const vector<double>& _data2, int window) {
    int ts_w = window + 1;
    rolling_data_container<> container(ts_w, 1);
    rolling_data_container<> container2(ts_w, 1);
    vector<double> row(1, 0);
    vector<double> row2(1, 0);
    ts_cross_rb_range nrdrr(1);
    nrdrr.set_row_size(window);

    for (int round = 0; round < 2; ++round) {
        container.clear();
        container2.clear();

        vector<double> x_, y;
        for (size_t i = 0; i < _data1.size(); ++i) {
            row[0] = _data1[i];
            container.push(row);
            row2[0] = _data2[i];
            container2.push(row2);
            nrdrr(container.get_old_row(), container2.get_old_row(), container.get_row_by_idx(ts_w - 1),
                  container2.get_row_by_idx(ts_w - 1), row.data());

            double _naive = dummy_ts_cross(_data1.data(), _data2.data(), i, window);
            REQUIRE(FloatEqual(row[0], _naive));
        }
    }
}

void test_ts_backward_cpn_rolling(const vector<double>& _data1, int window) {
    rolling_data_container<> container(window, 1);
    vector<double> row(1, 0);
    rolling_backward_cpn_rb_range nrdrr(1);
    nrdrr.set_row_size(window);

    for (int round = 0; round < 2; ++round) {
        nrdrr.init();
        container.clear();

        for (size_t i = 0; i < _data1.size(); ++i) {
            row[0] = _data1[i];
            container.push(row);
            double _naive = dummy_ts_backward_cpn(_data1.data(), i, window, 0);

            int real_window = window;
            if (container.m_count < window) {
                real_window = container.m_count;
            }

            {
                nrdrr.init();
                for (int ts_idx = 0; ts_idx < real_window; ++ts_idx) {
                    nrdrr.full_single(ts_idx, container.get_row_by_idx(ts_idx), row.data());
                }
                nrdrr.final_result(row.data());
                REQUIRE(FloatEqual(row[0], _naive));
            }
        }
    }
}

void test_ts_acp_rolling(const vector<double>& _data1, int window) {
    int lag = 1;
    rolling_data_container<> container(window, 1);
    vector<double> row(1, 0);
    rolling_ts_acp_rb_range nrdrr(1);
    nrdrr.set_row_size(window);

    for (int round = 0; round < 2; ++round) {
        nrdrr.init();
        container.clear();

        for (size_t i = 0; i < _data1.size(); ++i) {
            row[0] = _data1[i];
            container.push(row);
            double _naive = dummy_ts_acp(_data1.data(), i, window, lag);

            int real_window = window;
            if (container.m_count < window) {
                real_window = container.m_count;
            }

            {
                nrdrr.init();
                for (int ts_idx = 0; ts_idx < real_window - lag; ++ts_idx) {
                    nrdrr.full_single(ts_idx, container.get_row_by_idx(ts_idx), container.get_row_by_idx(ts_idx + lag));
                }
                nrdrr.final_result(row.data());
                REQUIRE(FloatEqual(row[0], _naive));
            }
        }
    }
}

void test_tsargmax_by_window(const vector<double>& _data1, int window) {
    int lag = 1;
    rolling_data_container<> container(window, 1);
    vector<double> row(1, 0);
    rolling_mq_percent_rb_range<double, std::less<double>> nrdrr(1);
    nrdrr.set_row_size(window);

    for (int round = 0; round < 2; ++round) {
        nrdrr.init();
        container.clear();

        for (size_t i = 0; i < _data1.size(); ++i) {
            row[0] = _data1[i];
            container.push(row);
            double _naive = dummy_ts_argmax(_data1.data(), i, window);
            nrdrr(container.get_old_row(), container.get_new_row(), row.data());
            if (!FloatEqual(row[0], _naive)) {
                _naive = dummy_ts_argmax(_data1.data(), i, window);
            }
            REQUIRE(FloatEqual(row[0], _naive));
        }
    }
}

void test_TsLteMean_method(int method, const vector<double>& data) {
    int window = 5, ins_num = 2;
    auto dummy_ret = ts_lte_mean(data, data, window, 0.5, NAN, method, 3, true);

    rolling_data_container<> container(window, ins_num);
    vector<double> row(ins_num, 0);
    TsLteMean rcs(ins_num);
    rcs.set_row_size(window);
    rcs.partial = true;
    rcs.method = method;

    for (int round = 0; round < 2; ++round) {
        container.clear();
        rcs.init();

        for (int j = 0; j < (int)data.size(); ++j) {
            double i = data[j];
            for (int k = 0; k < ins_num; ++k) {
                row[k] = i;
            }
            container.push(row);

            rcs(container.get_new_row(), container.get_new_row(), row.data());
            for (int k = 0; k < ins_num; ++k) {
                REQUIRE(FloatEqual(dummy_ret[j], row[k]));
            }
        }
    }
}

TEST_CASE("math_stats_rolling_rb_range", "[MathStatsRolling]") {
    vector<double> x = {1, 2, NAN, 3, NAN, NAN, 4, 5.5, 6, 7, NAN, 8, 9};
    vector<double> x1 = {3, 2, NAN, 1, NAN, NAN, 5, 3.5, 7, 6, NAN, 5, 7};
    vector<double> x2 = {2.5, 1.1, NAN, 0.7, NAN, NAN, 3.2, 1.4, 4.1, 3.5, NAN, 3.2, 4.2};
    vector<double> data{1, 2, 3, 4, 5, 6};
    vector<double> data2{6, 5, 4, 3, 2, 1};
    vector<double> weight{1, 1, 1, 1, 1, 1};
    uint32_t small_window = 4;
    uint32_t shift_count = data.size() - small_window;

    SECTION("variance online", "[MathStatsRolling]") {
        double var = 0, mean = 0;
        size_t count = 1;
        for (auto obv : data) ornate::variance_online(var, mean, obv, count);
        REQUIRE(var == 3.5);
        REQUIRE(mean == 3.5);
        REQUIRE(count == 7);

        ornate::variance_rolling<> vr(data.size());
        for (auto obv : data) vr(obv);
        REQUIRE(vr.variance == 3.5);
        REQUIRE(vr.mean == 3.5);
        REQUIRE(vr.count == 6);
    }

    SECTION("variance rolling", "[MathStatsRolling]") {
        vector<double> _data1{1, 2, 3, 4, 5, 6};
        test_variance_rolling(_data1, 3);
        test_variance_rolling(_data1, 5);

        vector<double> _data2{1, 2, NAN, 4, NAN, 6, 7, 8, 9, 10};
        test_variance_rolling(_data2, 3);
        test_variance_rolling(_data2, 5);

        vector<double> _data3{1.12, 1.13, NAN, 1.14, NAN, 1.12, 1.12, 1.13, 1.14, 1.12};
        test_variance_rolling(_data2, 3);
        test_variance_rolling(_data2, 5);
    }

    SECTION("covariance rolling", "[MathStatsRolling]") {
        test_covariance_rolling(data, data2, 3);
        test_covariance_rolling(data, data2, 5);

        vector<double> _data1{1, 2, NAN, 4, NAN, 6, 7, 8, 9, 10};
        vector<double> _data2{1.12, NAN, 1.13, 1.14, NAN, 1.12, 1.12, 1.13, 1.14, 1.12};
        test_covariance_rolling(_data1, _data2, 3);
        test_covariance_rolling(_data1, _data2, 5);
    }

    SECTION("covariance online", "[MathStatsRolling]") {
        double covariance = 0, mean1 = 0, mean2 = 0;
        size_t count = 1;
        for (size_t i = 0; i < data.size(); ++i) {
            ornate::covariance_online(covariance, mean1, mean2, data[i], data2[i], count);
        }
        REQUIRE(covariance == -3.5);
        REQUIRE(mean1 == 3.5);
        REQUIRE(mean2 == 3.5);
        REQUIRE(count == 7);

        ornate::covariance_rolling<> cvr(data.size());
        for (size_t i = 0; i < data.size(); ++i) {
            cvr(data[i], data2[i]);
        }
        REQUIRE(cvr.covariance == -3.5);
        REQUIRE(cvr.meanA == 3.5);
        REQUIRE(cvr.meanB == 3.5);
        REQUIRE(cvr.count == 6);
    }

    SECTION("covariance online small window", "[MathStatsRolling]") {
        double expected = -1.6666666666667;
        ornate::covariance_rolling<> cvr(small_window);
        for (size_t i = 0; i < data.size(); ++i) {
            cvr(data[i], data2[i]);
        }
        REQUIRE(FloatEqual(cvr.covariance, expected));
        REQUIRE(cvr.meanA == 4.5);
        REQUIRE(cvr.meanB == 2.5);
        REQUIRE(cvr.count == 4);

        REQUIRE(FloatEqual(covariance(data.data() + shift_count, data2.data() + shift_count, small_window), expected));
    }

    SECTION("corr online", "[MathStatsRolling]") {
        REQUIRE(FloatEqual(ornate::corr(data, data2), -1.0));
        REQUIRE(FloatEqual(ornate::cov(data, data2), -3.5));

        OnlineCorrelation oc;
        for (size_t i = 0; i < data.size(); ++i) {
            oc.Push(data[i], data2[i]);
        }
        REQUIRE(FloatEqual(oc.Result(), -1.0));
    }

    SECTION("corr rolling", "[MathStatsRolling]") {
        corr_rolling cr(data.size());

        double ret = 0;
        for (size_t i = 0; i < data.size(); ++i) {
            ret = cr(data[i], data2[i]);
        }
        REQUIRE(FloatEqual(ret, -1.0));

        test_corr_rolling(data, data2, 3);
        test_corr_rolling(data, data2, 5);

        vector<double> _data1{1, 2, NAN, 4, NAN, 6, 7, 8, 9, 10};
        vector<double> _data2{1.12, NAN, 1.13, 1.14, NAN, 1.12, 1.12, 1.13, 1.14, 1.12};
        test_corr_rolling(_data1, _data2, 3);
        test_corr_rolling(_data1, _data2, 5);
    }

    SECTION("skew rolling nan", "[MathStatsRolling]") {
        test_skew_by_window(x, 3);
        test_skew_by_window(x, 4);
        test_skew_by_window(x, 5);
        test_skew_by_window(x, 6);
    }

    SECTION("mean rolling", "[MathStatsRolling]") {
        int window = 4;
        rolling_mean_rb_range rmrr(1);
        rolling_mean_rb_range nrmrr(1);

        for (int round = 0; round < 2; ++round) {
            rmrr.init();
            nrmrr.init();

            mean_rolling<> sr(window);
            rolling_mean_rb rmr(window);
            rolling_data_container<> container(window, 1);
            vector<double> row(1, 0);
            double ret = 0;
            for (double i : data) {
                ret = sr(i);
                rmr(i);
                row[0] = i;
                container.push(row);
                rmrr(container.get_old_row(), container.get_new_row(), row.data());
                REQUIRE(FloatEqual(ret, rmr.mean));
                REQUIRE(FloatEqual(ret, row[0]));

                int real_window = window;
                if (container.m_count < window) {
                    real_window = container.m_count;
                }

                nrmrr.init();
                for (int ts_idx = 0; ts_idx < real_window; ++ts_idx) {
                    nrmrr.full_single(ts_idx, container.get_row_by_idx(ts_idx), row.data());
                }
                nrmrr.final_result(row.data());
                REQUIRE(FloatEqual(ret, row[0]));
            }
        }
    }

    SECTION("mean rolling nan", "[MathStatsRolling]") {
        int window = 4;
        rolling_data_container<> container(window, 1);
        vector<double> row(1, 0);
        rolling_mean_rb_range rmrr(1);
        rolling_mean_rb_range nrmrr(1);

        for (int round = 0; round < 2; ++round) {
            container.clear();
            rmrr.init();
            nrmrr.init();

            mean_rolling<> sr(window);
            rolling_mean_rb rmr(window);
            double ret = 0;
            for (double i : x) {
                ret = sr(i);
                rmr(i);
                row[0] = i;
                container.push(row);

                rmrr(container.get_old_row(), container.get_new_row(), row.data());
                REQUIRE(FloatEqual(ret, rmr.mean));
                REQUIRE(FloatEqual(ret, row[0]));

                int real_window = window;
                if (container.m_count < window) {
                    real_window = container.m_count;
                }

                nrmrr.init();
                for (int ts_idx = 0; ts_idx < real_window; ++ts_idx) {
                    nrmrr.full_single(ts_idx, container.get_row_by_idx(ts_idx), row.data());
                }
                nrmrr.final_result(row.data());
                REQUIRE(FloatEqual(ret, row[0]));
            }
        }
    }

    SECTION("kurtosis rolling nan", "[MathStatsRolling]") {
        test_kurtosis_by_window(x, 3);
        test_kurtosis_by_window(x, 4);
        test_kurtosis_by_window(x, 5);
        test_kurtosis_by_window(x, 6);
    }

    SECTION("decay rolling nan", "[MathStatsRolling]") {
        test_decay_by_window(x, 3);
        test_decay_by_window(x, 4);
        test_decay_by_window(x, 5);
        test_decay_by_window(x, 6);
    }

    SECTION("rank rolling nan", "[MathStatsRolling]") {
        vector<double> datum = {1, 4, NAN, NAN, 3, 1, 6, -2, 4, NAN, 7, 2, -3, NAN, NAN, 5};

        test_rank_by_window(datum, 3);
        test_rank_by_window(datum, 4);
        test_rank_by_window(datum, 5);
        test_rank_by_window(datum, 6);
    }

    SECTION("regression rolling nan", "[MathStatsRolling]") {
        test_regression_by_window(x, x1, 3);
        test_regression_by_window(x, x1, 4);
        test_regression_by_window(x, x1, 5);
        test_regression_by_window(x, x1, 6);
    }

    SECTION("regression3 rolling nan", "[MathStatsRolling]") {
        test_regression3_by_window(x, x1, x2, 3);
        test_regression3_by_window(x, x1, x2, 4);
        test_regression3_by_window(x, x1, x2, 5);
        test_regression3_by_window(x, x1, x2, 6);
    }

    SECTION("quantile rolling nan", "[MathStatsRolling]") {
        vector<double> datum = {1, 4, NAN, NAN, 3, 1, 6, -2, 4, NAN, 7, 2, -3, NAN, NAN, 5};

        test_quantile_by_window(datum, 3, 0.5);
        test_quantile_by_window(datum, 4, 0.5);
        test_quantile_by_window(datum, 5, 0.5);
        test_quantile_by_window(datum, 6, 0.5);
    }

    SECTION("ema_hl rolling nan", "[MathStatsRolling]") {
        vector<double> datum = {1, 4, NAN, NAN, 3, 1, 6, -2, 4, NAN, 7, 2, -3, NAN, NAN, 5};

        test_ema_hl_by_window(datum, 3);
        test_ema_hl_by_window(datum, 4);
        test_ema_hl_by_window(datum, 5);
        test_ema_hl_by_window(datum, 6);
    }

    SECTION("ema_hl_pp rolling nan", "[MathStatsRolling]") {
        vector<double> datum = {1, 4, NAN, NAN, 3, 1, 6, -2, 4, NAN, 7, 2, -3, NAN, NAN, 5};

        test_ema_hl_pp_by_window(datum, 3);
        test_ema_hl_pp_by_window(datum, 4);
        test_ema_hl_pp_by_window(datum, 5);
        test_ema_hl_pp_by_window(datum, 6);
    }

    SECTION("ols rolling nan", "[MathStatsRolling]") {
        test_ols_by_window(x, x1, 3);
        test_ols_by_window(x, x1, 4);
        test_ols_by_window(x, x1, 5);
        test_ols_by_window(x, x1, 6);
    }

    SECTION("ols3 rolling nan", "[MathStatsRolling]") {
        test_ols3_by_window(x, x1, x2, 3);
        test_ols3_by_window(x, x1, x2, 4);
        test_ols3_by_window(x, x1, x2, 5);
        test_ols3_by_window(x, x1, x2, 6);
    }

    SECTION("slope_no_intercept rolling nan", "[MathStatsRolling]") {
        test_slope_no_intercept_by_window(x1, 3);
        test_slope_no_intercept_by_window(x1, 4);
        test_slope_no_intercept_by_window(x1, 5);
        test_slope_no_intercept_by_window(x1, 6);
    }

    SECTION("slope rolling nan", "[MathStatsRolling]") {
        test_slope_by_window(x1, 3);
        test_slope_by_window(x1, 4);
        test_slope_by_window(x1, 5);
        test_slope_by_window(x1, 6);
    }

    SECTION("sharpe rolling nan", "[MathStatsRolling]") {
        test_sharpe_by_window(x1, 3);
        test_sharpe_by_window(x1, 4);
        test_sharpe_by_window(x1, 5);
        test_sharpe_by_window(x1, 6);
    }

    SECTION("scale rolling nan", "[MathStatsRolling]") {
        test_scale_by_window(x1, 3);
        test_scale_by_window(x1, 4);
        test_scale_by_window(x1, 5);
        test_scale_by_window(x1, 6);
    }

    SECTION("rank2 rolling nan", "[MathStatsRolling]") {
        test_rank2_by_window(x1, x2, 3);
        test_rank2_by_window(x1, x2, 4);
        test_rank2_by_window(x1, x2, 5);
        test_rank2_by_window(x1, x2, 6);
    }

    SECTION("dcor rolling", "[MathStatsRolling]") {
        test_dcor_rolling(data, data2, 3);
        test_dcor_rolling(data, data2, 5);

        vector<double> _data1{1, 2, NAN, 4, NAN, 6, 7, 8, 9, 10};
        vector<double> _data2{1.12, NAN, 1.13, 1.14, NAN, 1.12, 1.12, 1.13, 1.14, 1.12};
        test_dcor_rolling(_data1, _data2, 3);
        test_dcor_rolling(_data1, _data2, 5);
    }

    SECTION("tscross rolling", "[MathStatsRolling]") {
        vector<double> _data1{1, 2, NAN, 4, NAN, 1.11, 1.17, 1.11, 1.19, 1.05};
        vector<double> _data2{1.12, NAN, 1.13, 1.14, NAN, 1.12, 1.12, 1.13, 1.14, 1.12};
        test_tscross_rolling(_data1, _data2, 3);
        test_tscross_rolling(_data1, _data2, 5);
    }

    SECTION("ts_backward_cpn rolling", "[MathStatsRolling]") {
        vector<double> _data1{1, -2, NAN, 4, NAN, -1.11, 1.17, -1.11, -1.19, 1.05};
        test_ts_backward_cpn_rolling(_data1, 3);
        test_ts_backward_cpn_rolling(_data1, 4);
        test_ts_backward_cpn_rolling(_data1, 5);
    }

    SECTION("ts_acp rolling", "[MathStatsRolling]") {
        vector<double> _data1{1, -2, NAN, 4, NAN, -1.11, 1.17, -1.11, -1.19, 1.05};
        test_ts_acp_rolling(_data1, 3);
        test_ts_acp_rolling(_data1, 4);
        test_ts_acp_rolling(_data1, 5);
    }

    SECTION("tsargmax rolling nan", "[MathStatsRolling]") {
        vector<double> datum = {1, 4, NAN, NAN, 3, 1, 6, -2, 4, NAN, 7, 2, -3, NAN, NAN, 5};

        test_tsargmax_by_window(datum, 3);
        test_tsargmax_by_window(datum, 4);
        test_tsargmax_by_window(datum, 5);
        test_tsargmax_by_window(datum, 6);
    }

    SECTION("ts_gte_mean rolling nan", "[MathStatsRolling]") {
        int window = 5, ins_num = 2;
        auto dummy_ret = ts_gte_mean(data, data, window, 0.5, NAN, 1);
        std::vector<double> expected = {NAN, NAN, NAN, NAN, 4.0, 5.0};
        REQUIRE(FloatVecEqual(dummy_ret, expected));

        rolling_data_container<> container(window, ins_num);
        vector<double> row(ins_num, 0);
        TsGteMean rcs(ins_num);
        rcs.set_row_size(window);

        for (int round = 0; round < 2; ++round) {
            container.clear();
            rcs.init();

            for (int j = 0; j < (int)data.size(); ++j) {
                double i = data[j];
                for (int k = 0; k < ins_num; ++k) {
                    row[k] = i;
                }
                container.push(row);

                rcs(container.get_new_row(), container.get_new_row(), row.data());
                for (int k = 0; k < ins_num; ++k) {
                    REQUIRE(FloatEqual(expected[j], row[k]));
                }
            }
        }
    }

    SECTION("ts_lte_mean rolling nan", "[MathStatsRolling]") {
        int window = 5, ins_num = 2;
        auto dummy_ret = ts_lte_mean(data, data, window, 0.5, NAN, 1);
        std::vector<double> expected = {NAN, NAN, NAN, NAN, 2.0, 3.0};
        REQUIRE(FloatVecEqual(dummy_ret, expected));

        rolling_data_container<> container(window, ins_num);
        vector<double> row(ins_num, 0);
        TsLteMean rcs(ins_num);
        rcs.set_row_size(window);

        for (int round = 0; round < 2; ++round) {
            container.clear();
            rcs.init();

            for (int j = 0; j < (int)data.size(); ++j) {
                double i = data[j];
                for (int k = 0; k < ins_num; ++k) {
                    row[k] = i;
                }
                container.push(row);

                rcs(container.get_new_row(), container.get_new_row(), row.data());
                for (int k = 0; k < ins_num; ++k) {
                    REQUIRE(FloatEqual(expected[j], row[k]));
                }
            }
        }
    }

    SECTION("ts_lte_mean rolling nan partial method", "[MathStatsRolling]") {
        test_TsLteMean_method(0, data);
        test_TsLteMean_method(1, data);
        test_TsLteMean_method(2, data);
        test_TsLteMean_method(3, data);
        test_TsLteMean_method(4, data);
    }

    SECTION("ts_lte_mean rolling nan partial", "[MathStatsRolling]") {
        int window = 5, ins_num = 2;
        auto dummy_ret = ts_lte_mean(data, data, window, 0.5, NAN, 1, 3, true);
        std::vector<double> expected = {NAN, NAN, 1.5, 1.5, 2.0, 3.0};
        // REQUIRE(dummy_ret == expected);
        REQUIRE(FloatVecEqual(dummy_ret, expected));

        rolling_data_container<> container(window, ins_num);
        vector<double> row(ins_num, 0);
        TsLteMean rcs(ins_num);
        rcs.set_row_size(window);
        rcs.partial = true;

        for (int round = 0; round < 2; ++round) {
            container.clear();
            rcs.init();

            for (int j = 0; j < (int)data.size(); ++j) {
                double i = data[j];
                for (int k = 0; k < ins_num; ++k) {
                    row[k] = i;
                }
                container.push(row);

                rcs(container.get_new_row(), container.get_new_row(), row.data());
                for (int k = 0; k < ins_num; ++k) {
                    REQUIRE(FloatEqual(expected[j], row[k]));
                }
            }
        }
    }

    SECTION("ts_lte_max rolling nan partial", "[MathStatsRolling]") {
        int window = 5, ins_num = 2;
        auto dummy_ret = ts_lte_max(data, data, window, 0.5, NAN, 1, 3, true);
        std::vector<double> expected = {NAN, NAN, 2.0, 2.0, 3.0, 4.0};
        // REQUIRE(dummy_ret == expected);
        REQUIRE(FloatVecEqual(dummy_ret, expected));

        rolling_data_container<> container(window, ins_num);
        vector<double> row(ins_num, 0);
        TsLteMax rcs(ins_num);
        rcs.set_row_size(window);
        rcs.partial = true;

        for (int round = 0; round < 2; ++round) {
            container.clear();
            rcs.init();

            for (int j = 0; j < (int)data.size(); ++j) {
                double i = data[j];
                for (int k = 0; k < ins_num; ++k) {
                    row[k] = i;
                }
                container.push(row);

                rcs(container.get_new_row(), container.get_new_row(), row.data());
                for (int k = 0; k < ins_num; ++k) {
                    REQUIRE(FloatEqual(expected[j], row[k]));
                }
            }
        }
    }

    SECTION("ts_lte_min rolling nan partial", "[MathStatsRolling]") {
        int window = 5, ins_num = 2;
        auto dummy_ret = ts_lte_min(data, data, window, 0.5, NAN, 1, 3, true);
        std::vector<double> expected = {NAN, NAN, 1.0, 1.0, 1.0, 2.0};
        // REQUIRE(dummy_ret == expected);
        REQUIRE(FloatVecEqual(dummy_ret, expected));

        rolling_data_container<> container(window, ins_num);
        vector<double> row(ins_num, 0);
        TsLteMin rcs(ins_num);
        rcs.set_row_size(window);
        rcs.partial = true;

        for (int round = 0; round < 2; ++round) {
            container.clear();
            rcs.init();

            for (int j = 0; j < (int)data.size(); ++j) {
                double i = data[j];
                for (int k = 0; k < ins_num; ++k) {
                    row[k] = i;
                }
                container.push(row);

                rcs(container.get_new_row(), container.get_new_row(), row.data());
                for (int k = 0; k < ins_num; ++k) {
                    REQUIRE(FloatEqual(expected[j], row[k]));
                }
            }
        }
    }

    SECTION("ts_lte_sd rolling nan partial", "[MathStatsRolling]") {
        int window = 5, ins_num = 2;
        auto dummy_ret = ts_lte_sd(data, data, window, 0.5, NAN, 1, 3, true);
        std::vector<double> expected = {NAN, NAN, 0.7071067812, 0.7071067812, 1.0, 1.0};
        // REQUIRE(dummy_ret == expected);
        REQUIRE(FloatVecEqual(dummy_ret, expected));

        rolling_data_container<> container(window, ins_num);
        vector<double> row(ins_num, 0);
        TsLteSd rcs(ins_num);
        rcs.set_row_size(window);
        rcs.partial = true;

        for (int round = 0; round < 2; ++round) {
            container.clear();
            rcs.init();

            for (int j = 0; j < (int)data.size(); ++j) {
                double i = data[j];
                for (int k = 0; k < ins_num; ++k) {
                    row[k] = i;
                }
                container.push(row);

                rcs(container.get_new_row(), container.get_new_row(), row.data());
                for (int k = 0; k < ins_num; ++k) {
                    REQUIRE(FloatEqual(expected[j], row[k]));
                }
            }
        }
    }
}
