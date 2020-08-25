#include <math_stats_no_roll_rb_range.h>
#include <math_stats_rolling_rb_range.h>
#include <math_utils.h>
#include <iostream>
#include "catch.hpp"
#include "math_stats.h"

using namespace std;
using namespace ornate;

static vector<double> x = {1, 2, NAN, 3, NAN, NAN, 4, 5.5, 6, 7, NAN, 8, 9};
static vector<double> x1 = {3, 2, NAN, 1, NAN, NAN, 5, 3.5, 7, 6, NAN, 5, 7};
static vector<double> x2 = {2.5, 1.1, NAN, 0.7, NAN, NAN, 3.2, 1.4, 4.1, 3.5, NAN, 3.2, 4.2};
static vector<double> data{1, 2, 3, 4, 5, 6};
static vector<double> data2{6, 5, 4, 3, 2, 1};
static vector<double> weight{1, 1, 1, 1, 1, 1};
static uint32_t small_window = 4;
static uint32_t shift_count = data.size() - small_window;

TEST_CASE("variance online", "[MathStatsRolling]") {
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

void test_variance_rolling(const vector<double>& _data, int window) {
    ornate::variance_rolling<> vr(window);
    rolling_variance_rb rvr(window);
    rolling_data_container<> container(window, 1);
    vector<double> row(1, 0);
    rolling_variance_rb_range rvrr(1);
    no_roll_variance_rb_range nrvrr(1);

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
        if (std::isnan(std_naive)) {
            REQUIRE(std::isnan(s));
            REQUIRE(std::isnan(rvr.variance));
            REQUIRE(std::isnan(row[0]));
        } else {
            REQUIRE(sqrtf(vr.variance) == std_naive);
            REQUIRE(sqrtf(rvr.variance) == std_naive);
            REQUIRE(sqrtf(row[0]) == std_naive);
        }

        if (container.m_count >= window) {
            nrvrr.init();
            if (container.m_count >= window) {
                for (int ts_idx = 0; ts_idx < window; ++ts_idx) {
                    nrvrr(ts_idx, container.get_row_by_idx(ts_idx), row.data());
                }
                nrvrr.final_result(row.data());
            }
            if (std::isnan(std_naive)) {
                REQUIRE(std::isnan(row[0]));
            } else {
                REQUIRE(sqrtf(row[0]) == std_naive);
            }
        }
    }
}

TEST_CASE("variance rolling", "[MathStatsRolling]") {
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

void test_covariance_rolling(const vector<double>& _data1, const vector<double>& _data2, int window) {
    ornate::covariance_rolling<> cvr(window);
    rolling_cov_rb rcr(window);
    rolling_data_container<> container(window, 1);
    rolling_data_container<> container2(window, 1);
    vector<double> row(1, 0);
    vector<double> row2(1, 0);
    rolling_cov_rb_range rcrr(1);
    no_roll_cov_rb_range nrcrr(1);

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
        if (std::isnan(cov_naive)) {
            REQUIRE(std::isnan(s));
            REQUIRE(std::isnan(rcr.covariance));
            REQUIRE(std::isnan(row[0]));
        } else {
            REQUIRE(FloatEqual(cvr.covariance, cov_naive));
            REQUIRE(FloatEqual(rcr.covariance, cov_naive));
            REQUIRE(FloatEqual(row[0], cov_naive));
        }

        if (container.m_count >= window) {
            nrcrr.init();
            if (container.m_count >= window) {
                for (int ts_idx = 0; ts_idx < window; ++ts_idx) {
                    nrcrr(ts_idx, container.get_row_by_idx(ts_idx), container2.get_row_by_idx(ts_idx), row.data());
                }
                nrcrr.final_result(row.data());
            }
            if (std::isnan(cov_naive)) {
                REQUIRE(std::isnan(row[0]));
            } else {
                REQUIRE(FloatEqual(row[0], cov_naive));
            }
        }
    }
}

TEST_CASE("covariance rolling", "[MathStatsRolling]") {
    test_covariance_rolling(data, data2, 3);
    test_covariance_rolling(data, data2, 5);

    vector<double> _data1{1, 2, NAN, 4, NAN, 6, 7, 8, 9, 10};
    vector<double> _data2{1.12, NAN, 1.13, 1.14, NAN, 1.12, 1.12, 1.13, 1.14, 1.12};
    test_covariance_rolling(_data1, _data2, 3);
    test_covariance_rolling(_data1, _data2, 5);
}

TEST_CASE("covariance online", "[MathStatsRolling]") {
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

TEST_CASE("covariance online small window", "[MathStatsRolling]") {
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

TEST_CASE("corr online", "[MathStatsRolling]") {
    REQUIRE(FloatEqual(ornate::corr(data, data2), -1.0));
    REQUIRE(FloatEqual(ornate::cov(data, data2), -3.5));

    OnlineCorrelation oc;
    for (size_t i = 0; i < data.size(); ++i) {
        oc.Push(data[i], data2[i]);
    }
    REQUIRE(FloatEqual(oc.Result(), -1.0));
}

void test_corr_rolling(const vector<double>& _data1, const vector<double>& _data2, int window) {
    ornate::corr_rolling cvr(window);
    rolling_corr_rb rcr(window);
    rolling_data_container<> container(window, 1);
    rolling_data_container<> container2(window, 1);
    vector<double> row(1, 0);
    vector<double> row2(1, 0);
    rolling_corr_rb_range rcrr(1);
    no_roll_corr_rb_range nrcrr(1);

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

        if (std::isnan(corr_naive)) {
            REQUIRE(std::isnan(s));
            REQUIRE(std::isnan(rcr.corr));
            REQUIRE(std::isnan(row[0]));
        } else {
            REQUIRE(FloatEqual(s, corr_naive));
            REQUIRE(FloatEqual(rcr.corr, corr_naive));
            REQUIRE(FloatEqual(row[0], corr_naive));
        }

        if (container.m_count >= window) {
            nrcrr.init();
            if (container.m_count >= window) {
                for (int ts_idx = 0; ts_idx < window; ++ts_idx) {
                    nrcrr(ts_idx, container.get_row_by_idx(ts_idx), container2.get_row_by_idx(ts_idx), row.data());
                }
                nrcrr.final_result(row.data());
            }
            if (std::isnan(corr_naive)) {
                REQUIRE(std::isnan(row[0]));
            } else {
                REQUIRE(FloatEqual(row[0], corr_naive));
            }
        }
    }
}

TEST_CASE("corr rolling", "[MathStatsRolling]") {
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

static double calc_skew(const vector<double>& data_) {
    double mean_ = mean(data_);
    double std_ = 0;
    int valid_count = 0;
    for (double i : data_) {
        if (isfinite(i)) {
            ++valid_count;
            std_ += std::pow((i - mean_), 2);
        }
    }
    if (valid_count < 2) return NAN;
    std_ = sqrt(std_ / valid_count);
    if (std_ < 1e-7) return NAN;
    double ret = 0;
    for (double i : data_) {
        if (isfinite(i)) {
            ret += std::pow((i - mean_) / std_, 3);
        }
    }
    return ret / valid_count;
}

void test_skew_by_window(const vector<double>& x_, int window) {
    rolling_skew_rb sr(window);
    vector<double> y;

    rolling_data_container<> container(window, 1);
    vector<double> row(1, 0);
    rolling_skew_rb_range rsrr(1);
    no_roll_skew_rb_range nrsrr(1);

    double ret = 0;
    for (double d : x_) {
        if (y.size() < (size_t)window)
            y.push_back(d);
        else {
            for (int j = 1; j < window; ++j) {
                y[j - 1] = y[j];
            }
            y[window - 1] = d;
        }
        ret = sr(d);

        row[0] = d;
        container.push(row);
        rsrr(container.get_old_row(), container.get_new_row(), row.data());

        double expected = calc_skew(y);
        REQUIRE(FloatEqual(ret, expected));
        REQUIRE(FloatEqual(row[0], expected));

        if (container.m_count >= window) {
            nrsrr.init();
            if (container.m_count >= window) {
                for (int ts_idx = 0; ts_idx < window; ++ts_idx) {
                    nrsrr(ts_idx, container.get_row_by_idx(ts_idx), row.data());
                }
                nrsrr.final_result(row.data());
            }
            REQUIRE(FloatEqual(row[0], expected));
        }
    }
}

TEST_CASE("skew rolling nan", "[MathStatsRolling]") {
    test_skew_by_window(x, 3);
    test_skew_by_window(x, 4);
    test_skew_by_window(x, 5);
    test_skew_by_window(x, 6);
}

TEST_CASE("mean rolling", "[MathStatsRolling]") {
    int ts_window = 4;
    mean_rolling<> sr(ts_window);
    rolling_mean_rb rmr(ts_window);
    rolling_data_container<> container(ts_window, 1);
    vector<double> row(1, 0);
    rolling_mean_rb_range rmrr(1);
    no_roll_mean_rb_range nrmrr(1);

    double ret = 0;
    for (double i : data) {
        ret = sr(i);
        rmr(i);
        row[0] = i;
        container.push(row);
        rmrr(container.get_old_row(), container.get_new_row(), row.data());
        REQUIRE(FloatEqual(ret, rmr.mean));
        REQUIRE(FloatEqual(ret, row[0]));

        if (container.m_count >= ts_window) {
            nrmrr.init();
            if (container.m_count >= ts_window) {
                for (int ts_idx = 0; ts_idx < ts_window; ++ts_idx) {
                    nrmrr(ts_idx, container.get_row_by_idx(ts_idx), row.data());
                }
                nrmrr.final_result(row.data());
            }
            REQUIRE(FloatEqual(ret, row[0]));
        }
    }
}

TEST_CASE("mean rolling nan", "[MathStatsRolling]") {
    int ts_window = 4;
    mean_rolling<> sr(ts_window);
    rolling_mean_rb rmr(ts_window);
    rolling_data_container<> container(ts_window, 1);
    vector<double> row(1, 0);
    rolling_mean_rb_range rmrr(1);
    no_roll_mean_rb_range nrmrr(1);

    double ret = 0;
    for (double i : x) {
        ret = sr(i);
        rmr(i);
        row[0] = i;
        container.push(row);

        rmrr(container.get_old_row(), container.get_new_row(), row.data());
        REQUIRE(FloatEqual(ret, rmr.mean));
        REQUIRE(FloatEqual(ret, row[0]));

        if (container.m_count >= ts_window) {
            nrmrr.init();
            if (container.m_count >= ts_window) {
                for (int ts_idx = 0; ts_idx < ts_window; ++ts_idx) {
                    nrmrr(ts_idx, container.get_row_by_idx(ts_idx), row.data());
                }
                nrmrr.final_result(row.data());
            }
            REQUIRE(FloatEqual(ret, row[0]));
        }
    }
}

static double calc_kurtosis(const vector<double>& data_) {
    double mean_ = mean(data_);
    double std_ = 0;
    int valid_count = 0;
    for (double i : data_) {
        if (isfinite(i)) {
            ++valid_count;
            std_ += std::pow((i - mean_), 2);
        }
    }
    if (valid_count < 2) return NAN;
    std_ = sqrt(std_ / valid_count);
    if (std_ < 1e-7) return NAN;
    double ret = 0;
    for (double i : data_) {
        if (isfinite(i)) {
            ret += std::pow((i - mean_) / std_, 4);
        }
    }
    return ret / valid_count - 3.0;
}

void test_kurtosis_by_window(const vector<double>& x_, int window) {
    rolling_kurtosis_rb sr(window);
    vector<double> y;

    rolling_data_container<> container(window, 1);
    vector<double> row(1, 0);
    rolling_kurtosis_rb_range rkrr(1);
    no_roll_kurtosis_rb_range nrkrr(1);

    double ret = 0;
    for (double d : x_) {
        if (y.size() < (size_t)window)
            y.push_back(d);
        else {
            for (int j = 1; j < window; ++j) {
                y[j - 1] = y[j];
            }
            y[window - 1] = d;
        }
        ret = sr(d);

        row[0] = d;
        container.push(row);
        rkrr(container.get_old_row(), container.get_new_row(), row.data());

        double expected = calc_kurtosis(y);

        REQUIRE(FloatEqual(ret, expected));
        REQUIRE(FloatEqual(row[0], expected));

        if (container.m_count >= window) {
            nrkrr.init();
            if (container.m_count >= window) {
                for (int ts_idx = 0; ts_idx < window; ++ts_idx) {
                    nrkrr(ts_idx, container.get_row_by_idx(ts_idx), row.data());
                }
                nrkrr.final_result(row.data());
            }
            REQUIRE(FloatEqual(ret, row[0]));
        }
    }
}

TEST_CASE("kurtosis rolling nan", "[MathStatsRolling]") {
    test_kurtosis_by_window(x, 3);
    test_kurtosis_by_window(x, 4);
    test_kurtosis_by_window(x, 5);
    test_kurtosis_by_window(x, 6);
}

static double calc_decay(const vector<double>& data_) {
    int size = data_.size();
    double res = 0;
    int count = 0;
    for (int i = 0; i < size; ++i) {
        if (std::isfinite(data_[i])) {
            res += data_[i] * (i + 1);
            count += (i + 1);
        }
    }
    if (count > 0) {
        return res / count;
    } else {
        return NAN;
    }
}

void test_decay_by_window(const vector<double>& x_, int window) {
    vector<double> y;

    rolling_data_container<> container(window, 1);
    vector<double> row(1, 0);
    rolling_decay_rb_range rdrr(1);
    rdrr.set_row_size(window);
    no_roll_decay_rb_range nrdrr(1);
    nrdrr.set_row_size(window);

    for (double d : x_) {
        if (y.size() < (size_t)window)
            y.push_back(d);
        else {
            for (int j = 1; j < window; ++j) {
                y[j - 1] = y[j];
            }
            y[window - 1] = d;
        }

        row[0] = d;
        container.push(row);
        rdrr(container.get_old_row(), container.get_new_row(), row.data());

        double expected = calc_decay(y);
        if (!FloatEqual(row[0], expected)) {
            cout << row[0] << " " << expected << endl;
            rdrr(container.get_old_row(), container.get_new_row(), row.data());
            cout << row[0] << " " << expected << endl;
        }
        REQUIRE(FloatEqual(row[0], expected));

        if (container.m_count >= window) {
            nrdrr.init();
            if (container.m_count >= window) {
                for (int ts_idx = 0; ts_idx < window; ++ts_idx) {
                    nrdrr(ts_idx, container.get_row_by_idx(ts_idx), row.data());
                }
                nrdrr.final_result(row.data());
            }
            REQUIRE(FloatEqual(expected, row[0]));
        }
    }
}

TEST_CASE("decay rolling nan", "[MathStatsRolling]") {
    test_decay_by_window(x, 3);
    test_decay_by_window(x, 4);
    test_decay_by_window(x, 5);
    test_decay_by_window(x, 6);
}

void test_rank_by_window(const vector<double>& x_, int window) {
    rolling_rank_rb<double> rrrb(window);
    vector<double> y;

    rolling_data_container<> container(window, 2);
    vector<double> row(2, 0);
    rolling_rank_rb_range<double> rrrr(2);
    rrrr.set_row_size(window);

    double ret = 0;
    for (double d : x_) {
        if (y.size() < (size_t)window)
            y.push_back(d);
        else {
            for (int j = 1; j < window; ++j) {
                y[j - 1] = y[j];
            }
            y[window - 1] = d;
        }
        ret = rrrb(d);

        row[0] = d;
        row[1] = d;
        container.push(row);
        rrrr(container.get_old_row(), container.get_new_row(), row.data());

        vector<double> tmp = y;
        ornate::rank(tmp);
        double expected = tmp.back();

        REQUIRE(FloatEqual(ret, expected));
        REQUIRE(FloatEqual(row[0], expected));
        REQUIRE(FloatEqual(row[1], expected));
    }
}

TEST_CASE("rank rolling nan", "[MathStatsRolling]") {
    vector<double> datum = {1, 4, NAN, NAN, 3, 1, 6, -2, 4, NAN, 7, 2, -3, NAN, NAN, 5};

    test_rank_by_window(datum, 3);
    test_rank_by_window(datum, 4);
    test_rank_by_window(datum, 5);
    test_rank_by_window(datum, 6);
}

void test_regression_by_window(const vector<double>& x_, const vector<double>& y_, int window) {
    regression2_rolling rr(window);
    vector<double> _x, _y;
    rolling_regression2_rb rrrb(window);

    rolling_data_container<> container(window, 2);
    rolling_data_container<> container2(window, 2);
    vector<double> row(2, 0);
    vector<double> row2(2, 0);
    rolling_regression2_rb_range rrrr(2);
    rrrr.set_row_size(window);
    no_roll_regression2_rb_range nrrrr(2);
    nrrrr.set_row_size(window);

    double ret = 0;
    for (size_t i = 0; i < x_.size(); ++i) {
        if (_x.size() < (size_t)window) {
            _x.push_back(x_[i]);
            _y.push_back(y_[i]);
        } else {
            for (int j = 1; j < window; ++j) {
                _x[j - 1] = _x[j];
                _y[j - 1] = _y[j];
            }
            _x[window - 1] = x_[i];
            _y[window - 1] = y_[i];
        }

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
        rrrr(container.get_old_row(), container2.get_old_row(), container.get_new_row(), container2.get_new_row(),
             row.data(), row2.data());
        REQUIRE(FloatEqual(row[0], a));
        REQUIRE(FloatEqual(row[1], a));
        REQUIRE(FloatEqual(row2[0], b));
        REQUIRE(FloatEqual(row2[1], b));

        if (container.m_count >= window) {
            nrrrr.init();
            if (container.m_count >= window) {
                for (int ts_idx = 0; ts_idx < window; ++ts_idx) {
                    nrrrr(ts_idx, container.get_row_by_idx(ts_idx), container2.get_row_by_idx(ts_idx));
                }
                nrrrr.final_result(row.data(), row2.data());
            }
            REQUIRE(FloatEqual(a, row[0]));
            REQUIRE(FloatEqual(b, row2[0]));
        }
    }
}

TEST_CASE("regression rolling nan", "[MathStatsRolling]") {
    test_regression_by_window(x, x1, 3);
    test_regression_by_window(x, x1, 4);
    test_regression_by_window(x, x1, 5);
    test_regression_by_window(x, x1, 6);
}

void test_regression3_by_window(const vector<double>& x1_, const vector<double>& x2_, const vector<double>& y_,
                                int window) {
    vector<double> _x1, _x2, _y;
    rolling_regression3_rb rrrb(window);

    rolling_data_container<> container(window, 2);
    rolling_data_container<> container2(window, 2);
    rolling_data_container<> container3(window, 2);
    vector<double> row(2, 0);
    vector<double> row2(2, 0);
    vector<double> row3(2, 0);
    rolling_regression3_rb_range rrrr(2);
    rrrr.set_row_size(window);
    no_roll_regression3_rb_range nrrrr(2);
    nrrrr.set_row_size(window);

    double ret = 0;
    for (size_t i = 0; i < x1_.size(); ++i) {
        if (_x1.size() < (size_t)window) {
            _x1.push_back(x1_[i]);
            _x2.push_back(x2_[i]);
            _y.push_back(y_[i]);
        } else {
            for (int j = 1; j < window; ++j) {
                _x1[j - 1] = _x1[j];
                _x2[j - 1] = _x2[j];
                _y[j - 1] = _y[j];
            }
            _x1[window - 1] = x1_[i];
            _x2[window - 1] = x2_[i];
            _y[window - 1] = y_[i];
        }

        double b0, b1, b2;
        ornate::regression3(_y, _x1, _x2, &b0, &b1, &b2);

        rrrb(y_[i], x1_[i], x2_[i]);

        if (!FloatEqual(rrrb.b0, b0)) {
            cout << rrrb.b0 << endl;
        }
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
             container2.get_new_row(), container3.get_new_row(), row.data(), row2.data(), row3.data());
        REQUIRE(FloatEqual(row[0], b0));
        REQUIRE(FloatEqual(row[1], b0));
        REQUIRE(FloatEqual(row2[0], b1));
        REQUIRE(FloatEqual(row2[1], b1));
        REQUIRE(FloatEqual(row3[0], b2));
        REQUIRE(FloatEqual(row3[1], b2));

        if (container.m_count >= window) {
            nrrrr.init();
            if (container.m_count >= window) {
                for (int ts_idx = 0; ts_idx < window; ++ts_idx) {
                    nrrrr(ts_idx, container.get_row_by_idx(ts_idx), container2.get_row_by_idx(ts_idx),
                          container3.get_row_by_idx(ts_idx));
                }
                nrrrr.final_result(row.data(), row2.data(), row3.data());
            }
            REQUIRE(FloatEqual(b0, row[0]));
            REQUIRE(FloatEqual(b1, row2[0]));
            REQUIRE(FloatEqual(b2, row3[0]));
        }
    }
}

TEST_CASE("regression3 rolling nan", "[MathStatsRolling]") {
    test_regression3_by_window(x, x1, x2, 3);
    test_regression3_by_window(x, x1, x2, 4);
    test_regression3_by_window(x, x1, x2, 5);
    test_regression3_by_window(x, x1, x2, 6);
}

void test_quantile_by_window(const vector<double>& x_, int window, double percent) {
    rolling_quantile_rb<double> rrrb(window, percent);
    vector<double> y;

    rolling_data_container<> container(window, 2);
    vector<double> row(2, 0);
    rolling_quantile_rb_range<double> rqrr(2, percent);
    rqrr.set_row_size(window);

    double ret = 0;
    for (double d : x_) {
        if (y.size() < (size_t)window)
            y.push_back(d);
        else {
            for (int j = 1; j < window; ++j) {
                y[j - 1] = y[j];
            }
            y[window - 1] = d;
        }
        ret = rrrb(d);

        row[0] = d;
        row[1] = d;
        container.push(row);
        rqrr(container.get_old_row(), container.get_new_row(), row.data());

        vector<double> tmp = y;
        double expected = ornate::quantile(tmp, percent);

        if (!FloatEqual(ret, expected)) {
            cout << ret << " " << expected << endl;
        }
        REQUIRE(FloatEqual(ret, expected));
        REQUIRE(FloatEqual(row[0], expected));
        REQUIRE(FloatEqual(row[1], expected));
    }
}

TEST_CASE("quantile rolling nan", "[MathStatsRolling]") {
    vector<double> datum = {1, 4, NAN, NAN, 3, 1, 6, -2, 4, NAN, 7, 2, -3, NAN, NAN, 5};

    test_quantile_by_window(datum, 3, 0.5);
    test_quantile_by_window(datum, 4, 0.5);
    test_quantile_by_window(datum, 5, 0.5);
    test_quantile_by_window(datum, 6, 0.5);
}

void test_ema_hl_by_window(const vector<double>& x_, int window) {
    rolling_ema_hl_rb rrrb(window, 4);
    vector<double> y;

    rolling_data_container<> container(window, 2);
    vector<double> row(2, 0);
//    rolling_quantile_rb_range<double> rqrr(2, percent);
//    rqrr.set_row_size(window);

    no_roll_ema_hl_rb_range nrqrr(2);
    nrqrr.set_row_size(window);
    nrqrr.set_param("hl", "4");

    double ret = 0;
    for (double d : x_) {
        if (y.size() < (size_t)window)
            y.push_back(d);
        else {
            for (int j = 1; j < window; ++j) {
                y[j - 1] = y[j];
            }
            y[window - 1] = d;
        }
        ret = rrrb(d);

//        row[0] = d;
//        row[1] = d;
//        container.push(row);
//        rqrr(container.get_old_row(), container.get_new_row(), row.data());

        vector<double> tmp = y;
        double expected = ornate::ema_hl(tmp, window, window - 1, 4);

        if (!FloatEqual(ret, expected)) {
            cout << ret << " " << expected << endl;
        }
        REQUIRE(FloatEqual(ret, expected));
//        REQUIRE(FloatEqual(row[0], expected));
//        REQUIRE(FloatEqual(row[1], expected));

        if (container.m_count >= window) {
            nrqrr.init();
            if (container.m_count >= window) {
                for (int ts_idx = 0; ts_idx < window; ++ts_idx) {
                    nrqrr(ts_idx, container.get_row_by_idx(ts_idx), row.data());
                }
                nrqrr.final_result(row.data());
            }
            REQUIRE(FloatEqual(expected, row[0]));
        }
    }
}

TEST_CASE("ema_hl rolling nan", "[MathStatsRolling]") {
    vector<double> datum = {1, 4, NAN, NAN, 3, 1, 6, -2, 4, NAN, 7, 2, -3, NAN, NAN, 5};

    test_ema_hl_by_window(datum, 3);
    test_ema_hl_by_window(datum, 4);
    test_ema_hl_by_window(datum, 5);
    test_ema_hl_by_window(datum, 6);
}
