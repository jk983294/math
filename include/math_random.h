#ifndef ORNATE_MATH_RANDOM_H
#define ORNATE_MATH_RANDOM_H

#include <random>
#include "math_utils.h"

using std::mt19937;
using std::normal_distribution;
using std::random_device;
using std::uniform_int_distribution;
using std::uniform_real_distribution;
using std::vector;

namespace ornate {

vector<uint32_t> generate_stream_data(size_t count) {
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
}  // namespace ornate

#endif
