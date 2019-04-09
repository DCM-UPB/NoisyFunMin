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
// Luckily, especially if one restricts the class to support solely the push_back operation, i.e. no
// pop() or arbitrary insert/remove, the implementation comes down to a few methods consisting of simple
// single if-clauses (with the exception of the reserve()/setcap() methods).
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
template <class ValueT> // in NFM, ValueT is usually NoisyValue, NoisyIOPair or just int/double(vectors)
class PushBackBuffer
{
protected:
    std::vector<ValueT> _vec; // the vector for actually storing the elements
    size_t _ncap; // the buffer capacity (max number of elements)
    size_t _inext; // the next index to be written at on push

    void _incINext() noexcept { _inext = (++_inext < _ncap) ? _inext : _inext - _ncap; } // this is a bit cheaper than modulo
    size_t _wrapIndex(size_t i) const noexcept { return (i < _ncap) ? i : i - _ncap; }

public:
    explicit PushBackBuffer(size_t size = 0/*initial full capacity*/) noexcept;

    // --- Element access (const)

    const ValueT &front() const; // oldest element
    const ValueT &back() const; // newest element
    const ValueT &operator[](size_t i) const; // insertion-time-based access (i.e. [0] == front, [size()] == back)

    // Raw access to underlying data
    // Use these whenever you need to read the data without
    // requiring information about relative position in time!
    const std::vector<ValueT> & vec() const noexcept { return _vec; } // e.g. averaging elements from vec().begin() to vec().end() is always correct
    const ValueT * data() const noexcept { return _vec.data(); } // similarly, going from data() to data()+size() would be correct

    // --- Capacity

    bool empty() const noexcept { return _vec.empty(); };
    bool full() const noexcept { return (_vec.size() == _ncap); }

    size_t size() const noexcept { return _vec.size(); } // this is max(npushed, ncap)
    size_t capacity() const noexcept { return _ncap; }

    void reserve(size_t new_cap); // resembles the usual meaning of reserve() (i.e. may only increase cap)
    void setcap(size_t new_cap); // also allows to shrink (preserving newest elements)

    // --- push_back and emplace_back, clear and swap

    void push_back(const ValueT &val); // copy version
    void push_back(ValueT &&val); // move version

    template <class... Args>
    void emplace_back(Args &&... args); // construct in-place version

    void clear() noexcept; // reduce size to 0, leave capacity as is
    void swap(PushBackBuffer &other) noexcept;
};

// non-member swap
template<class ValueT>
void swap( PushBackBuffer<ValueT>& lhs, PushBackBuffer<ValueT>& rhs ) noexcept { lhs.swap(rhs); }


// --- Implementation

template <class ValueT>
PushBackBuffer<ValueT>::PushBackBuffer(const size_t size) noexcept {
    _vec.reserve(size); // we pre-allocate the required space
    _ncap = size;
    _inext = 0;
}

template <class ValueT>
const ValueT &PushBackBuffer<ValueT>::front() const {
    return (_inext < _vec.size()) ? _vec[_inext] : _vec.front(); // second case only if empty (UB) or not full
}

template <class ValueT>
const ValueT &PushBackBuffer<ValueT>::back() const {
    return (_inext != 0) ? _vec[_inext - 1] : _vec.back(); // works always (UB when empty)
}

template <class ValueT>
const ValueT &PushBackBuffer<ValueT>::operator[](const size_t i) const { // when index i out of bounds, UB
    return (_inext < _vec.size()) ? _vec[this->_wrapIndex(_inext + i)] : _vec[i];
}

template <class ValueT>
void PushBackBuffer<ValueT>::reserve(const size_t new_cap)
{
    if (new_cap <= _ncap) { return; } // like for std::vector, reserve does nothing if new_cap<_ncap

    if (_inext < _vec.size()) {
        std::vector<ValueT> new_vec; // the easiest way seems to be using a new vector
        new_vec.reserve(new_cap);
        for (size_t i = _inext; i < _vec.size(); ++i) { // fill new vec, starting from oldest
            new_vec.emplace_back(std::move(_vec[i])); // we can move elements
        }
        for (size_t i = 0; i < _inext; ++i) { // now fill in the rest
            new_vec.emplace_back(std::move(_vec[i]));
        }
        _vec.swap(new_vec);
        _inext = _vec.size(); // there is at least one unwritten element at the end of vec
    }
    else { // _inext == _vec.size();
        _vec.reserve(new_cap);
    }

    _ncap = new_cap;
}

template <class ValueT>
void PushBackBuffer<ValueT>::setcap(const size_t new_cap)
{
    if (new_cap>=_ncap) { // simply reuse reserve (helps to reduce cases)
        this->reserve(new_cap);
        return;
    }

    if (_inext < _vec.size()) { // vec is full and will be overfull after cap reduce
        const size_t dcap = _ncap - new_cap; // number of elements to remove, > 0
        std::vector<ValueT> new_vec;
        new_vec.reserve(new_cap);
        for (size_t i = _inext + dcap; i < _vec.size(); ++i) {
            new_vec.emplace_back(std::move(_vec[i])); // we can move elements
        }
        const size_t drest = (_inext + dcap < _vec.size()) ? 0 : _inext + dcap + 1 - _vec.size(); // rest to skip
        for (size_t i = drest; i < _inext; ++i) { // now fill in the rest
            new_vec.emplace_back(std::move(_vec[i]));
        }
        _vec.swap(new_vec);
    }
    else if (_vec.size() > new_cap) { // not full before but overfull after cap reduce
        const size_t dsize = _vec.size() - new_cap; // number of elements to remove, > 0
        std::vector<ValueT>(_vec.begin()+dsize, _vec.end()).swap(_vec); // remove the first dsize elements
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
    this->_incINext();
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
    this->_incINext();
}

template <class ValueT>
template <class... Args>
void PushBackBuffer<ValueT>::emplace_back(Args &&... args) {
    if (this->full()) {
        _vec[_inext] = ValueT(std::forward<Args>(args)...);
    }
    else {
        _vec.emplace_back(std::forward<Args>(args)...);
    }
    this->_incINext();
}

template <class ValueT>
void PushBackBuffer<ValueT>::clear() noexcept
{
    _vec.clear();
    _inext = 0;
}

template <class ValueT>
void PushBackBuffer<ValueT>::swap(PushBackBuffer &other) noexcept {
    _vec.swap(other._vec);
    std::swap(_ncap, other._ncap);
    std::swap(_inext, other._inext);
}
} // namespace nfm

#endif
