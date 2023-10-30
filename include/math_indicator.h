#pragma once


#include <cmath>
#include <cstddef>

namespace ornate {
/**
 * from Kenneth H. Shaleen - Volume And Open Interest
 */
double HHS(const double* close, const double* volume, const double* oi, std::size_t i, std::size_t n);
}