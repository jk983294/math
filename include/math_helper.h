#ifndef ORNATE_MATH_HELPER_H
#define ORNATE_MATH_HELPER_H

namespace ornate {
/**
 * suffix dt means different type of arguments which differ from std::plus, it only support same type
 */
template <typename T1, typename T2, typename TRet = T1>
struct plus_dt {
    TRet operator()(const T1& x, const T2& y) const { return x + y; }
};

template <typename T1, typename T2, typename TRet = T1>
struct minus_dt {
    TRet operator()(const T1& x, const T2& y) const { return x - y; }
};

template <typename T1, typename T2, typename TRet = T1>
struct divide_dt {
    TRet operator()(const T1& x, const T2& y) const { return x / y; }
};

template <typename T1, typename T2, typename TRet = T1>
struct multiply_dt {
    TRet operator()(const T1& x, const T2& y) const { return x * y; }
};
}  // namespace ornate

#endif
