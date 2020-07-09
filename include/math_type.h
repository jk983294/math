#ifndef __MATH_TYPE_H_
#define __MATH_TYPE_H_

#include <cmath>
#include <cstddef>
#include <cstdint>

#define IN
#define OUT
#define INOUT

namespace ornate {

typedef int32_t int32;
typedef uint32_t uint32;
typedef char char8;
typedef int16_t int16;
typedef uint16_t uint16;
typedef char int8;
typedef uint8_t uint8;
typedef uint8_t uchar8;
typedef int64_t int64;
typedef uint64_t uint64;
typedef float float32;
typedef double double64;

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

bool inline isvalid(float n) { return std::isfinite(n); }
bool inline isvalid(double n) { return std::isfinite(n); }
bool inline isvalid(uint32_t n) { return (n != uint32_nan); }
bool inline isvalid(int32_t n) { return (n != int32_nan); }
bool inline isvalid(uint16_t n) { return (n != uint16_nan); }
bool inline isvalid(int16_t n) { return (n != int16_nan); }
bool inline isvalid(uint8_t n) { return (n != uint8_nan); }
bool inline isvalid(int8_t n) { return (n != int8_nan); }
bool inline isvalid(bool n) { return ((uint8)n != bool_nan); }
bool inline isvalid(const char *n) { return (n != nullptr && n[0] != 0); }

template <typename T>
inline constexpr T get_nan() {
    return 0;
}
template <>
inline constexpr float get_nan<float>() {
    return NAN;
}
template <>
inline constexpr double get_nan<double>() {
    return NAN;
}
template <>
inline constexpr int32_t get_nan<int32_t>() {
    return int32_nan;
}
template <>
inline constexpr uint32_t get_nan<uint32_t>() {
    return uint32_nan;
}

}  // namespace ornate

#endif
