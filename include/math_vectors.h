#ifndef ORNATE_MATH_VECTORS_H
#define ORNATE_MATH_VECTORS_H

#include "math_utils.h"

using std::vector;

namespace ornate {
/**
 * T should be some integer type
 */
template <typename T>
inline double kendallsTau(const vector<T>& rankX, const vector<T>& rankY) {
    if (rankX.size() != rankY.size()) {
        throw std::runtime_error("Input vectors must have the same size.");
    }

    int n = rankX.size();
    int concordant = 0, discordant = 0;

    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 1; j < n; ++j) {
            int xDiff = rankX[i] - rankX[j];
            int yDiff = rankY[i] - rankY[j];

            if (xDiff * yDiff > 0) {
                concordant++; // Concordant pair
            } else if (xDiff * yDiff < 0) {
                discordant++; // Discordant pair
            }
            // Ignore ties (xDiff == 0 or yDiff == 0)
        }
    }

    double totalPairs = n * (n - 1) / 2.0;
    double tau = (concordant - discordant) / totalPairs;
    return tau;
}
}  // namespace ornate

#endif
