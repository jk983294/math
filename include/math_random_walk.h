#ifndef ORNATE_MATH_RANDOM_WALK_H
#define ORNATE_MATH_RANDOM_WALK_H

#include <random>
#include "math_utils.h"

using std::max;
using std::min;

namespace ornate {
/**
 * general brownian motion
 * @param price_path
 * @param start_price
 * @param N         number of time steps
 * @param r         riskless interest rate
 * @param q         divident yield
 * @param sigma     volatility of stock
 * @param T         time (expiry)
 */
void GBM(std::vector<double>& price_path, double start_price, const int N = 1000, const double r = 0.001,
         const double q = 0.0, const double sigma = 0.00010, const double T = 1) {
    std::random_device rd;
    std::mt19937 generator(rd());
    std::normal_distribution<double> dist(0.0, 1.0);
    double dt = T / N;
    price_path.clear();
    double drift = exp(dt * ((r - q) - 0.5 * sigma * sigma));
    double vol = sqrt(sigma * sigma * dt);
    for (int i = 0; i < N; i++) {
        double Z = dist(generator);
        start_price = start_price * drift * exp(vol * Z);
        price_path.push_back(start_price);
    }
}

/**
 * generate open, high, low from close_prices
 */
void GBM_ohlc(std::vector<double>& close_prices, std::vector<double>& open, std::vector<double>& high,
              std::vector<double>& low, double up_limit_percent = 1.1, double low_limit_percent = 0.9) {
    std::random_device rd;
    std::mt19937 generator(rd());
    size_t len = close_prices.size();

    open.resize(len);
    high.resize(len);
    low.resize(len);

    double range = up_limit_percent - low_limit_percent;
    std::normal_distribution<double> up_gen(low_limit_percent + 0.6 * range, range / 4);
    std::normal_distribution<double> down_gen(low_limit_percent + 0.4 * range, range / 4);
    std::normal_distribution<double> open_gen(low_limit_percent + 0.5 * range, range / 2);

    for (size_t i = 0; i < len; i++) {
        double open_ = open_gen(generator);
        double up_ = up_gen(generator);
        double down_ = down_gen(generator);
        high[i] = min(max(max(up_, down_), open_), up_limit_percent) * close_prices[i];
        low[i] = max(min(min(up_, down_), open_), low_limit_percent) * close_prices[i];
        open_ = min(open_, up_limit_percent);
        open_ = max(open_, low_limit_percent);
        open[i] = open_ * close_prices[i];
    }
}
}  // namespace ornate

#endif
