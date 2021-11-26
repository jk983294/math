#ifndef ORNATE_MATH_RING_H
#define ORNATE_MATH_RING_H

#include <math_vector.h>
#include "math_utils.h"

namespace ornate {

template <class T>
class RingArray {
public:
    explicit RingArray(std::size_t n)
        : _n(n),
          _capacity(n < 2 ? 2 : ((n & (n - 1)) == 0 ? n : roundup_pow_of_two(n))),
          _capacity_sub_one(_capacity - 1) {
        _data.resize(_capacity, T());
    }

    ~RingArray() = default;

    bool empty() const { return _size == 0; }

    bool full() const { return _size == _n; }

    // Assume ring is not full, otherwise UB
    void push_back(const T& val) {
        ++_total;
        ++_size;
        if (_size > 1) {
            _tail = (_tail + 1) & _capacity_sub_one;
        } else {
            _head = _tail;
        }
        _data[_tail] = val;
    }

    // Assume ring is not empty, otherwise UB
    void pop_front() {
        --_size;
        _head = (_head + 1) & _capacity_sub_one;
    }

    // Assume ring has at least 2 elements, otherwise UB
    void pop_front_push_back(const T& val) {
        ++_total;
        _head = (_head + 1) & _capacity_sub_one;
        _tail = (_tail + 1) & _capacity_sub_one;
        _data[_tail] = val;
    }

    // Assume ring has at least 2 elements, otherwise UB
    T pushOverride(const T& val) {
        T old = _data[_tail];
        _data[_tail] = val;
        return old;
    }

    // Assume ring is not empty, otherwise UB
    void pop_back() {
        --_size;
        _tail = (_tail - 1) & _capacity_sub_one;
    }

    // safe
    void push(const T& val) {
        if (full()) {
            pop_front_push_back(val);
        } else {
            push_back(val);
        }
    }

    const T& operator[](std::size_t cbar) const { return _data[(cbar - _total + _tail) & _capacity_sub_one]; }

    const T& get_by_offset(int offset) const { return _data[(_tail - offset) & _capacity_sub_one]; }

    const T& front() const { return _data[_head]; }

    const T& back() const { return _data[_tail]; }

    std::size_t size() const { return _size; }

    std::size_t capacity() const { return _capacity; }

    int total() const { return _total; }
    int nbar() const { return _total + 1; }

    void clear() {
        _total = -1;
        _size = _head = _tail = 0;
    }

private:
    std::size_t roundup_pow_of_two(std::size_t n) {
        n |= (n >> 1);
        n |= (n >> 2);
        n |= (n >> 4);
        n |= (n >> 8);
        n |= (n >> 16);
        return n + 1;
    }

private:
    std::size_t _n{0};                 // window size
    std::size_t _capacity{0};          // max size
    std::size_t _capacity_sub_one{0};  // _capacity - 1
    int _total = -1;
    std::size_t _size{};  // actual size
    std::size_t _head{};
    std::size_t _tail{};
    std::vector<T> _data;
};

}  // namespace ornate

#endif
