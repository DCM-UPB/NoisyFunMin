#ifndef NFM_PUSHBACKBUFFER_HPP
#define NFM_PUSHBACKBUFFER_HPP

#include <vector>

namespace nfm
{

template <class ValueT>
class PushBackBuffer
{
protected:
    std::vector<ValueT> _vec; // the vector for actually storing the elements
    size_t _ncap; // the buffer capacity (max number of elements)
    size_t _inext; // the next index to be written at on push

    size_t _incIndex(size_t i) const noexcept { return (++i < _ncap) ? i : i - _ncap; }
    size_t _wrapIndex(size_t i) const noexcept { return (i < _ncap) ? i : i - _ncap; }
    void _incINext() noexcept { _inext = (++_inext < _ncap) ? _inext : _inext - _ncap; }

public:
    explicit PushBackBuffer(size_t size) noexcept
    {
        _vec.reserve(size); // we pre-allocate the required space
        _ncap = size;
        _inext = 0;
    }


    // --- Element access (const)

    const ValueT &front() const
    {
        return (_inext < _vec.size()) ? _vec[_inext] : _vec.front(); // second case only if empty (UB) or not full
    }

    const ValueT &back() const
    {
        return (_inext != 0) ? _vec[_inext - 1] : _vec.back(); // works always (UB when empty)
    }

    const ValueT &operator[](size_t i) const
    {
        return (_inext < _vec.size()) ? _vec[this->_wrapIndex(_inext + i)] : _vec[i];
    }

    // --- Capacity

    bool empty() const noexcept { return _vec.empty(); };
    bool full() const noexcept { return (_vec.size() == _ncap); }

    size_t size() const noexcept { return _vec.size(); }
    size_t max_size() const noexcept { return _vec.max_size(); }
    size_t capacity() const noexcept { return _ncap; }

    void reserve(size_t new_cap) // somewhat resembles the usual meaning of reserve()
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
            _vec = std::move(new_vec);
        }
        else {
            _vec.reserve(new_cap);
        }
        _inext = _ncap;
        _ncap = new_cap;
    }

    /*
    void setCapacity(size_t new_cap) // also allows to shrink (preserving newest elements)
    {
        if (new_cap>=_ncap) { this->reserve(new_cap); }
        const size_t diff = _ncap - new_cap;

        std::vector<ValueT> new_vec;
        new_vec.reserve(new_cap);

        if (_inext < _vec.size()) {
            for (size_t i = _inext + diff; i < _vec.size(); ++i) {
                new_vec.emplace_back(std::move(_vec[i])); // we can move elements
            }
            const size_t drest = (_inext + diff < _vec.size()) ? 0 : _inext + diff + 1 - _vec.size(); // rest to skip
            for (size_t i = drest; i < _inext; ++i) { // now fill in the rest
                new_vec.emplace_back(std::move(_vec[i]));
            }
        }
        else {
            for (size_t i = diff; i< _vec.size(); ++i) {
                new_vec.emplace_back(std::move(_vec[i]));
            }
        }

        _vec = std::move(new_vec);
        _inext =

    }*/

    // --- push_back and emplace_back, clear and swap

    void push_back(const ValueT &val) // push ref
    {
        if (this->full()) {
            _vec[_inext] = val; // use copy assignment
        }
        else {
            _vec.push_back(val);
        }
        this->_incINext();
    }

    void push_back(ValueT &&val) // push rvref
    {
        if (this->full()) {
            _vec[_inext] = val; // use move assignment
        }
        else {
            _vec.push_back(val);
        }
        this->_incINext();
    }

    template <class... Args>
    void emplace_back(Args &&... args)
    {
        if (this->full()) {
            _vec[_inext] = ValueT(std::forward<Args>(args)...);
        }
        else {
            _vec.emplace_back(std::forward<Args>(args)...);
        }
        this->_incINext();
    }

    void clear() noexcept // capacity remains!
    {
        _vec.clear();
        _inext = 0;
    }

    void swap( PushBackBuffer & other )
    {
        _vec.swap(other);
        std::swap(_ncap, other._ncap);
        std::swap(_inext, other._inext);
    }
};
} // namespace nfm

#endif
