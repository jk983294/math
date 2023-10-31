#include <math_lm.h>
#include <math_random.h>
#include <math_stats.h>
#include <math_vector.h>
#include <armadillo>
#include <cmath>
#include <cstdio>
#include <stats.hpp>
#include <tuple>

namespace ornate {

void Model::fit_lm(bool intercept) {
    m_lm.total_row = total_row;
    m_lm.m_features = m_features;
    m_lm.m_y = m_y;
    m_lm.m_intercept = intercept;

    m_lm.lm(false);
}

}  // namespace ornate