#pragma once

#include <cmath>
#include <cstddef>

namespace ornate {
/**
 * from Kenneth H. Shaleen - Volume And Open Interest
 * seems good for financial instrument, not good for commodity
 */
double HHS(const double* close, const double* volume, const double* oi, std::size_t i, std::size_t n,
           double blowup_factor = 1.5);
}  // namespace ornate