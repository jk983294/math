#ifndef __MATH_TYPE_H_
#define __MATH_TYPE_H_

#include <cmath>
#include <cstdint>

namespace ornate {

constexpr int32_t int32_nan = 0x80000000;
constexpr uint32_t uint32_nan = 0xffffffff;
constexpr uint16_t uint16_nan = 0xffff;
constexpr int16_t int16_nan = 0x8000;
constexpr uint8_t bool_nan = 0xff;
constexpr uint8_t uint8_nan = 0xff;
constexpr int8_t int8_nan = 0x80;
constexpr float float32_nan = NAN;
constexpr double double64_nan = NAN;
constexpr double float64_nan = NAN;

}  // namespace ornate

#endif
