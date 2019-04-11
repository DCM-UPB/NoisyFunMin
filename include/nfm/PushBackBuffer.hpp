#ifndef NFM_PUSHBACKBUFFER_HPP
#define NFM_PUSHBACKBUFFER_HPP

#include <vector>

namespace nfm
{

// Simplified (push-only/auto-pop) circular buffer
// Author: Jan Kessler (2019)
//
// In many numeric applications that follow an iterative scheme, one may have the desire to keep
// track of the results of a certain amount of previous iterations, without permanently storing
// every single step of the whole program loop. This means that on every step the new calculated
// value should be stored and the oldest stored value should be dropped.
// Surprisingly, none of the standard STL containers allows to do this to full satisfaction. Using
// a plain std::vector allows fast iteration over the elements, but the pop_front operation leads
// to a lot of copy operations per step. Using a std::list avoids the copying, but iterating over
// it is ugly and slow (+ memory is all over the place). std::deque would be the still unsatisfactory
// compromise of both.
// But, in this special use case of a fixed-size buffer, one may use an overall superior data structure
// instead: A ring buffer. This buffer stores the values in contiguous memory internally, but allows
// insertion-order-based indexing in O(1) complexity with very small overhead. Pushing a new element
// comes down to a single copy or move operation. Unfortunately, ring buffer implementations are only
// available via external libraries like boost and we don't want to needlessly add such dependency.
// But especially if one restricts the class to support solely the push_back operation, i.e. no pop()
// or arbitrary insert/remove, the implementation comes down to a few methods consisting of simple
// single if-clauses (with the exception of the reserve()/set_cap() methods). Furthermore, we only
// need to provide const versions of accessor methods. New values enter the buffer only via push/emplace.
//
// NOTE: The class is not a fully STL-compatible container (i.e. fulfills not all requirements), but
// the public interface largely resembles a subset of the std::vector interface. Quite notably though,
// we don't provide iterator methods (e.g. begin(), end()), restricting the direct use of the container
// inside STL algorithms. However, whenever the actual position of the buffer elements do not matter,
// you may simply use the underlying vector directly (via buffer.vec()) or even the memory itself (buffer.data()).
//
// NOTE 2: In most of my other code I don't care about noexcept specifiers, but the std::vector interface
//         uses it, so for conformity I chose to use noexcept consistently here.
//
//
// Benchmark: http://quick-bench.com/DrWoTCBHvwJ7u8jdbpC2LB-2Q20
//
template <class ValueT> // in NFM, ValueT is usually NoisyValue, NoisyIOPair or just int/double(vectors)
class PushBackBuffer
{
protected:
    std::vector<ValueT> _vec; // the vector for actually storing the elements
    size_t _ncap; // the buffer capacity (max number of elements)
    size_t _inext; // the next internal storage index to be written at on push

    void _inc_inext() noexcept; // increment and wrap inext

public:
    // Constructor (NOTE: When default-constructed, any accessor or push is UB)
    explicit PushBackBuffer(size_t size = 0/*initial full capacity*/) noexcept;

    // --- Element access (const)

    const ValueT &front() const; // oldest element
    const ValueT &back() const; // newest element
    const ValueT &operator[](size_t i) const; // insertion-time-based access (i.e. [0] == front, [size()] == back)

    // Raw access to underlying data
    // Use these whenever you need to read the data without
    // requiring information about relative position in time!
    const std::vector<ValueT> &vec() const noexcept { return _vec; } // e.g. averaging elements from vec().begin() to vec().end() is always correct
    const ValueT * data() const noexcept { return _vec.data(); } // similarly, going from data() to data()+size() would be correct

    // --- Capacity

    bool empty() const noexcept { return _vec.empty(); }
    bool full() const noexcept { return (_vec.size() == _ncap); } // i.e. cap=0 -> is full

    size_t size() const noexcept { return _vec.size(); } // this is max(npushed, ncap)
    size_t capacity() const noexcept { return _ncap; }

    void reserve(size_t new_cap); // resembles the usual meaning of reserve() (i.e. may only increase cap)
    void set_cap(size_t new_cap); // also allows to shrink (preserving newest elements)

    // --- push_back and emplace_back, clear and swap

    void push_back(const ValueT &val); // copy version
    void push_back(ValueT &&val); // move version

    template <class... Args>
    void emplace_back(Args &&... args); // construct in-place version

    void clear() noexcept; // reduce size to 0, leave capacity as is
    void swap(PushBackBuffer &other) noexcept;
};

// non-member swap
template <class ValueT>
void swap(PushBackBuffer<ValueT> &lhs, PushBackBuffer<ValueT> &rhs) noexcept { lhs.swap(rhs); }


// --- Implementation

template <class ValueT>
PushBackBuffer<ValueT>::PushBackBuffer(const size_t size) noexcept
{
    _vec.reserve(size); // we pre-allocate the required space
    _ncap = size;
    _inext = 0;
}

template <class ValueT>
void PushBackBuffer<ValueT>::_inc_inext() noexcept
{
    if (++_inext == _ncap) { _inext = 0; } // this beats modulo and is safe here
}

template <class ValueT>
const ValueT &PushBackBuffer<ValueT>::front() const
{
    return (_inext < _vec.size()) ? _vec[_inext] : _vec.front(); // second case only if empty (UB) or not full
}

template <class ValueT>
const ValueT &PushBackBuffer<ValueT>::back() const
{
    return (_inext != 0) ? _vec[_inext - 1] : _vec.back(); // works always (UB when empty)
}

template <class ValueT>
const ValueT &PushBackBuffer<ValueT>::operator[](const size_t i) const
{ // when index i out of bounds, UB
    return (_inext < _vec.size()) ? _vec[(_inext + i)%_ncap] : _vec[i];
}

template <class ValueT>
void PushBackBuffer<ValueT>::reserve(const size_t new_cap)
{
    if (new_cap <= _ncap) { return; } // like for std::vector, reserve does nothing if new_cap<_ncap

    if (_inext < _vec.size()) {
        std::vector<ValueT> new_vec; // the easiest way seems to be using a new vector
        new_vec.reserve(new_cap);
        for (size_t i = 0; i < _vec.size(); ++i) {
            new_vec.push_back(std::move(_vec[(_inext + i)%_ncap])); // we can move the elements
        }
        _vec.swap(new_vec);
        _inext = _vec.size(); // there is space for at least one new element at the end of vec
    }
    else { // _inext == _vec.size();
        _vec.reserve(new_cap);
    }

    _ncap = new_cap;
}

template <class ValueT>
void PushBackBuffer<ValueT>::set_cap(size_t new_cap)
{
    if (new_cap >= _ncap) { // simply reuse reserve (helps to reduce cases)
        this->reserve(new_cap);
        return;
    }

    if (_inext < _vec.size()) { // vec is full and will be overfull after cap reduce
        const size_t dcap = _ncap - new_cap; // number of elements to remove, > 0
        std::vector<ValueT> new_vec;
        new_vec.reserve(new_cap);
        for (size_t i = dcap; i < _vec.size(); ++i) {
            new_vec.push_back(std::move(_vec[(_inext + i)%_ncap])); // we can move the elements
        }
        _vec.swap(new_vec);
    }
    else if (_vec.size() > new_cap) { // not full before but overfull after cap reduce
        const size_t dsize = _vec.size() - new_cap; // number of elements to remove, > 0
        std::vector<ValueT>(_vec.begin() + dsize, _vec.end()).swap(_vec); // remove the first dsize elements
    }
    else { // not full before and _vec.size() <= new_cap
        _vec.shrink_to_fit(); // else the vector cap will not decrease
        _vec.reserve(new_cap);
    }

    _ncap = new_cap;
    if (this->full()) { _inext = 0; }  // at this point _vec is always ordered, so 0 is oldest
}

template <class ValueT>
void PushBackBuffer<ValueT>::push_back(const ValueT &val)
{
    if (this->full()) {
        _vec[_inext] = val; // use copy assignment
    }
    else {
        _vec.push_back(val);
    }
    this->_inc_inext();
}

template <class ValueT>
void PushBackBuffer<ValueT>::push_back(ValueT &&val)
{
    if (this->full()) {
        _vec[_inext] = val; // use move assignment
    }
    else {
        _vec.push_back(val);
    }
    this->_inc_inext();
}

template <class ValueT>
template <class... Args>
void PushBackBuffer<ValueT>::emplace_back(Args &&... args)
{
    if (this->full()) {
        _vec[_inext] = ValueT(std::forward<Args>(args)...); // I'm not sure what exactly this compiles to
    }
    else {
        _vec.emplace_back(std::forward<Args>(args)...);
    }
    this->_inc_inext();
}

template <class ValueT>
void PushBackBuffer<ValueT>::clear() noexcept
{
    _vec.clear();
    _inext = 0;
}

template <class ValueT>
void PushBackBuffer<ValueT>::swap(PushBackBuffer &other) noexcept
{
    _vec.swap(other._vec);
    std::swap(_ncap, other._ncap);
    std::swap(_inext, other._inext);
}
} // namespace nfm

#endif
