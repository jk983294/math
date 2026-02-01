#pragma once

#include <vector>
#include <algorithm>
#include "math_type.h"
#include "math_vector.h"

namespace ornate {
/**
 * fill daily data(val) into intra data(out)
 * ukey2ii is d_univ's mapping, if not provided, built first, then use it to fill
 */
inline void daily2intra(const int32_t* d_univ, const int32_t* i_univ,
  const double* val, double* out, int d_len, int i_len,
  const std::unordered_map<int, int>* ukey2ii = nullptr, double fill = NAN) {
  std::fill(out, out + i_len, fill);
  
  std::unordered_map<int, int> daily_ukey_to_ii;
  if (ukey2ii == nullptr) {
    for (int ii = 0; ii < d_len; ii++) {
      int daily_ukey = d_univ[ii];
      daily_ukey_to_ii[daily_ukey] = ii;
    }
    ukey2ii = &daily_ukey_to_ii;
  }
  
  for (int i_idx = 0; i_idx < i_len; i_idx++) {
    int intra_ukey = i_univ[i_idx];
    auto itr = ukey2ii->find(intra_ukey);
    if (itr != ukey2ii->end()) {
      out[i_idx] = val[itr->second];
    }
  }
}
}