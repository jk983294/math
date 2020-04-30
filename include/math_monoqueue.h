#include <math_utils.h>
#include <functional>

namespace ornate {

/**
 * default is min queue
 */
template<typename TData, typename TCmp=std::greater<TData> >
struct MonotoneQueue {
    struct Cell {
        int seq;
        TData data;
    };
    
    const int capacity = 0;
    int front = 0;
    int rear = 0;
    int seq = -1;
    Cell* _data = nullptr;
    TCmp cmp;

    MonotoneQueue(int size): capacity{size + 1}, cmp() {
        _data = new Cell[capacity];
    }

    ~MonotoneQueue() {
        delete[] _data;
    }

    Cell* operator[](int i){
        return _data + i;
    }

    Cell& GetCell(int i) {
        return _data[i];
    }

    const Cell& GetCell(int i) const {
        return _data[i];
    }

    const TData& GetCellData(int ind){
        return _data[ind].data;
    }

    int GetCellSeq(int ind){
        return _data[ind].seq;
    }

    void Push(const TData& value) {
        // pop outdated data
        int oldest_seq = seq - capacity + 2; // slide window left boundary
        int ptr = front;
        for (; ptr != rear && GetCellSeq(ptr) <= oldest_seq; ptr = (ptr + 1) % capacity);
        front = ptr; // take the front to the left boundary, (lazy delete)
        
        // push new data
        ptr = (rear + capacity - 1) % capacity;
        int end_ptr = (front + capacity - 1) % capacity;
        for (; ptr != end_ptr && cmp(GetCellData(ptr), value); ptr = (ptr + capacity - 1) % capacity);
        ptr = (ptr + 1) % capacity; // push new data, and lazy-delete (move right boundary leftwards) larger items

        auto& cell = GetCell(ptr);
        cell.seq = ++seq;
        cell.data = value;
        rear = (ptr + 1) % capacity;
    }

    TData Top() {
        if(empty()) return get_nan<TData>();
        return GetCell(front).data;
    }

    int TopIndex() {
        if(empty()) return -1;
        return seq - GetCell(front).seq;
    }

    bool empty() { return front == rear; }

    void Clear() {
        front = 0;
        rear = 0;
        seq = -1;
    }
};

}
