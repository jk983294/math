#ifndef ORNATE_MATH_RANDOM_H
#define ORNATE_MATH_RANDOM_H

#include <algorithm>
#include <random>
#include "math_utils.h"

using std::mt19937;
using std::normal_distribution;
using std::random_device;
using std::uniform_int_distribution;
using std::uniform_real_distribution;
using std::vector;

namespace ornate {

inline vector<uint32_t> generate_stream_data(size_t count) {
    random_device rd;
    mt19937 generator(rd());
    std::uniform_int_distribution<uint32_t> uid(0, 1000000);
    vector<uint32_t> ret;
    ret.resize(count);
    for (size_t i = 0; i < count; ++i) {
        ret[i] = uid(generator);
    }
    return ret;
}

template <typename T, template <typename> class Dist>
vector<T> generate_random_data(size_t count, T a, T b) {
    random_device rd;  // non-deterministic generator
    mt19937 generator(rd());
    Dist<T> uid(a, b);
    vector<T> ret;
    ret.resize(count);
    for (size_t i = 0; i < count; ++i) {
        ret[i] = uid(generator);
    }
    return ret;
}

template <typename T, template <typename> class Dist>
vector<vector<T>> generate_random_data(size_t m, size_t n, T a, T b) {
    random_device rd;  // non-deterministic generator
    mt19937 generator(rd());
    Dist<T> uid(a, b);
    vector<vector<T>> ret;
    ret.resize(m);
    for (size_t i = 0; i < m; ++i) {
        ret[i].resize(n);
        for (size_t j = 0; j < n; ++j) {
            ret[i][j] = uid(generator);
        }
    }
    return ret;
}

template <typename T>
vector<T> generate_uniform_int(size_t count, T from = 0, T to = 1000000) {
    return generate_random_data<T, uniform_int_distribution>(count, from, to);
}

template <typename T>
vector<T> generate_uniform_float(size_t count, T from = 0.0, T to = 1000000.0) {
    return generate_random_data<T, uniform_real_distribution>(count, from, to);
}

template <typename T>
vector<T> generate_gaussian(size_t count, T mean = 0.0, T delta = 1.0) {
    return generate_random_data<T, normal_distribution>(count, mean, delta);
}

template <typename T = double>
inline vector<vector<T>> generate_gaussian_matrix(size_t m, size_t n, T mean = 0, T delta = 1) {
    return generate_random_data<T, normal_distribution>(m, n, mean, delta);
}

template <typename T = double>
inline vector<vector<T>> generate_uniform_matrix(size_t m, size_t n, T from = 0, T to = 1) {
    return generate_random_data<T, uniform_real_distribution>(m, n, from, to);
}

template <typename T = double>
inline vector<double> build_choice_vector(const vector<T> &weight) {
    vector<double> ret;
    double sum_ = std::accumulate(weight.begin(), weight.end(), 0.0);
    double accum = 0;
    for (auto val : weight) {
        accum += val;
        ret.push_back(accum / sum_);
    }
    return ret;
}

inline int choice_with_accum_weight(const vector<double> &accum_weight, double prob) {
    auto itr = std::upper_bound(accum_weight.begin(), accum_weight.end(), prob);
    if (itr == accum_weight.end())
        return (int)accum_weight.size() - 1;
    else
        return itr - accum_weight.begin();
}

// select k from n
inline std::vector<int> reservoir_sampling(int n, int k) {
    std::vector<int> reservoir(k);
    std::iota(reservoir.begin(), reservoir.end(), 0);
    if (k >= n) {
        reservoir.erase(reservoir.begin() + n, reservoir.end());
    } else {
        mt19937 generator(std::random_device{}());
        std::uniform_int_distribution<int> toReplace(0, static_cast<int>(k));
        for (int i = k; i < n; i++) {
            std::bernoulli_distribution whether2replace{static_cast<double>(k) / i};
            if (whether2replace(generator)) reservoir[toReplace(generator)] = i;
        }
        std::sort(reservoir.begin(), reservoir.end());
    }
    return reservoir;
}

inline std::vector<int> sample_by_n(int total, int n) { return reservoir_sampling(total, n); }

inline std::vector<int> sample_by_ratio(int total, double ratio) {
    return reservoir_sampling(total, static_cast<int>(total * ratio));
}

struct MyRandom {
    MyRandom() : generator(std::random_device()()), urd(0., 1.0) {}

    MyRandom(int from, int to) : generator(std::random_device()()), uid(from, to) {}

    double random() { return urd(generator); }

    int random_int(int from, int to) {
        std::uniform_int_distribution<int> uid_(from, to);
        return uid_(generator);
    }

    int random_int() { return uid(generator); }

    template <typename T>
    T random_from(const std::vector<T> &datum) {
        if (datum.empty())
            return T();
        else if (datum.size() == 1)
            return datum.front();
        else {
            int idx = random_int(0, (int)datum.size() - 1);
            return datum[idx];
        }
    }

    template <typename T>
    T choice(const std::vector<T> &datum, const std::vector<double> &accum_weight) {
        if (datum.empty())
            return T();
        else if (datum.size() == 1)
            return datum.front();
        else {
            double prob = random();
            int idx = choice_with_accum_weight(accum_weight, prob);
            return datum[idx];
        }
    }

    vector<int> choice_no_replace(int len, int n) {
        std::vector<int> ret(len, 0);
        std::iota(ret.begin(), ret.end(), 0);
        if (n < len) {
            std::shuffle(ret.begin(), ret.end(), generator);
            return vector<int>(ret.begin(), ret.begin() + n);
        } else {
            return ret;
        }
    }

    std::mt19937 generator;
    std::uniform_real_distribution<double> urd;
    std::uniform_int_distribution<int> uid;
};
}  // namespace ornate

#endif
