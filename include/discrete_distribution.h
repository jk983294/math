#ifndef ORNATE_DISCRETE_DISTRIBUTION_H_
#define ORNATE_DISCRETE_DISTRIBUTION_H_

#include <algorithm>
#include <cstddef>
#include <initializer_list>
#include <random>
#include <unordered_set>
#include <vector>

namespace ornate {
class discrete_weights {
public:
    using pointer = double const*;
    using iterator = double const*;

    discrete_weights() = default;

    /**
     * weights = Weight values. The weights must be non-negative finite numbers.
     * Construct a complete binary tree in which the leaves contain the weights
     * and internal nodes contain the sums of the weights of children.
     */
    explicit discrete_weights(std::vector<double> const& weights) {
        std::size_t leaves = 1;  // Determine the number of leaves.

        for (;;) {
            if (leaves >= weights.size()) {
                break;
            }
            leaves *= 2;
        }

        // We store the binary tree as an array using the usual scheme:
        //
        //   root = 0 ,
        //   parent(node) = (node - 1) / 2 .
        _sumtree.resize(leaves * 2 - 1);
        _leaves = leaves;
        _events = weights.size();

        // Fill the leaves.
        for (std::size_t i = 0; i < _events; i++) {
            if (std::isfinite(weights[i])) {
                _sumtree[_leaves + i - 1] = weights[i];
                ++_effective_n;
            } else
                _sumtree[_leaves + i - 1] = 0;
        }

        // Then, fill internal nodes from leaves to the root
        for (std::size_t layer = 1;; layer++) {
            auto const layer_size = leaves >> layer;
            if (layer_size == 0) {
                break;
            }

            auto const start = layer_size - 1;
            auto const end = start + layer_size;

            for (std::size_t node = start; node < end; node++) {
                auto const lchild = 2 * node + 1;
                auto const rchild = 2 * node + 2;
                _sumtree[node] = _sumtree[lchild] + _sumtree[rchild];
            }
        }
    }

    discrete_weights(std::initializer_list<double> const& weights) : discrete_weights{std::vector<double>{weights}} {}

    /**
     * Returns the number of events.
     */
    [[nodiscard]] std::size_t size() const noexcept { return _events; }

    /**
     * Returns a pointer to the array containing weight values.
     */
    [[nodiscard]] pointer data() const noexcept { return _sumtree.data() + _leaves - 1; }

    /**
     * Returns an iterator pointing to the beginning of the array containing weight values.
     */
    [[nodiscard]] iterator begin() const noexcept { return data(); }

    /**
     * Returns an iterator pointing to the past the end of the array containing weight values.
     */
    [[nodiscard]] iterator end() const noexcept { return data() + size(); }

    /**
     * Returns the weight of the i-th event.
     */
    inline double operator[](std::size_t i) const { return *(data() + i); }

    /**
     * Returns the sum of the weights.
     */
    [[nodiscard]] double sum() const { return _sumtree[0]; }

    /**
     * Updates the weight of the i-th event.
     *
     * Behavior is undefined if `i` is out of range or `weight` is
     * negative or not finite. It is also undefined that the sum of
     * weights overflow due to the update.
     *
     * Params:
     *   i      = Index of the event to update weight.
     *   weight = New weight value.
     */
    void update(std::size_t i, double weight) {
        auto node = _leaves + i - 1;
        _sumtree[node] = weight;

        do {
            node = (node - 1) / 2;
            auto const lchild = 2 * node + 1;
            auto const rchild = 2 * node + 2;
            _sumtree[node] = _sumtree[lchild] + _sumtree[rchild];
        } while (node > 0);
    }

    /**
     * Finds the event whose cumulative weight interval covers given probe value.
     *
     * Let `w[i]` be the weight of the i-th event, and let `s[i]` be the
     * cumulative sum until the i-th event. That is,
     *
     *    s[0] = 0 ,
     *    s[i] = w[0] + ... + w[i-1] .
     *
     * Then, this function finds the event i such that
     *
     *    s[i] <= probe < s[i+1] .
     *
     * Params:
     *   probe = Probe weight used to find an event.
     *
     * Returns:
     *   The index of the event found.
     */
    [[nodiscard]] std::size_t find(double probe) const {
        std::size_t node = 0;

        for (;;) {
            auto const lchild = 2 * node + 1;
            auto const rchild = 2 * node + 2;

            if (lchild >= _sumtree.size()) {
                break;
            }

            if (probe < _sumtree[lchild]) {
                node = lchild;
            } else {
                probe -= _sumtree[lchild];
                node = rchild;
            }
        }

        auto index = node - (_leaves - 1);

        // Search may overshoot due to numerical errors.
        if (index >= _events) {
            index = _events - 1;
        }

        return index;
    }

    std::size_t effective_n() const { return _effective_n; }

private:
    std::vector<double> _sumtree;
    std::size_t _leaves = 0;
    std::size_t _events = 0;
    std::size_t _effective_n = 0;
};

inline bool operator==(ornate::discrete_weights const& w1, ornate::discrete_weights const& w2) {
    if (w1.size() != w2.size()) {
        return false;
    }
    return std::equal(w1.begin(), w1.end(), w2.begin());
}

inline bool operator!=(ornate::discrete_weights const& w1, ornate::discrete_weights const& w2) { return !(w1 == w2); }

/**
 * Distribution of random integers with given weights.
 */
class discrete_distribution {
public:
    /**
     * Class holding distribution parameter set, namely, the weights.
     */
    class param_type : public ornate::discrete_weights {
    public:
        using ornate::discrete_weights::discrete_weights;
        param_type(ornate::discrete_weights const& weights) : discrete_weights{weights} {}
    };

    discrete_distribution() = default;

    /**
     * Creates a discrete distribution over `[0, N)` with given weights where `N = weights.size()`.
     *
     * Params:
     *   weights = Weight values. The weights must be non-negative finite numbers.
     */
    explicit discrete_distribution(std::vector<double> const& weights) : _weights{weights} {}
    discrete_distribution(std::initializer_list<double> const& weights) : _weights{weights} {}
    explicit discrete_distribution(param_type const& param) : _weights{param} {}
    param_type const& param() const noexcept { return _weights; }

    /**
     * Reconfigures the distribution using given parameter set
     */
    void param(param_type const& new_param) { _weights = new_param; }

    /**
     * Returns the minimum possible integer generated from this
     * distribution, namely, zero.
     */
    size_t min() const { return 0; }

    /**
     * Returns the maximum possible integer generated from this distribution.
     */
    size_t max() const { return _weights.size() - 1; }

    /**
     * Returns the sum of the weights.
     */
    double sum() const { return _weights.sum(); }

    /**
     * Updates the weight of the number `i`.
     * Params:
     *   i      = Number to change weight. Must be in the valid interval `[min(), max()]`.
     *   weight = New weight. Must be non-negative finite number.
     */
    void update(size_t i, double weight) { return _weights.update(i, weight); }

    /**
     * Generates an integer randomly from the weighted distribution.
     *
     * Params:
     *   random = Random number generator to use.
     *
     * Time complexity:
     *   O(log N) where N is the upper bound.
     */
    template <typename RNG>
    size_t operator()(RNG& random) const {
        std::uniform_real_distribution<> dis(0.0, _weights.sum());
        double val = dis(random);
        return _weights.find(val);
    }

    std::size_t effective_n() const { return _weights.effective_n(); }

private:
    param_type _weights;
};

inline std::vector<std::size_t> weighted_sample(const std::vector<double>& weight, std::size_t n,
                                                bool replace = false) {
    std::vector<std::size_t> ret;
    ret.reserve(n);
    ornate::discrete_distribution distr(weight);
    if (!replace && n >= distr.effective_n()) {
        ret.resize(weight.size());
        std::iota(ret.begin(), ret.end(), 0);
    } else {
        std::random_device rd;
        std::mt19937 gen(rd());

        if (replace) {
            while (ret.size() < n) {
                ret.push_back(distr(gen));
            }
        } else {
            std::unordered_set<std::size_t> s;
            while (s.size() < n) {
                size_t choice = distr(gen);
                if (s.find(choice) == s.end()) {
                    s.insert(choice);
                }
            }
            std::copy(s.begin(), s.end(), std::back_inserter(ret));
        }
    }
    return ret;
}
}  // namespace ornate

#endif
