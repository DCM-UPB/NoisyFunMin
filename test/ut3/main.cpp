#include <cassert>
#include <iostream>
#include <numeric>

#include "nfm/PushBackBuffer.hpp"
#include "nfm/LogManager.hpp"

#include "TestNFMFunctions.hpp"

bool verbose = false; // to enable cout output

template <class ValueT>
void assertBuffer(const nfm::PushBackBuffer<ValueT> &buffer, size_t ncap_expect, size_t nsize_expect, const ValueT * range, const ValueT * data)
{
    if (verbose) { std::cout << "nsize " << buffer.size() << ", nsize_expect " << nsize_expect << std::endl; }
    if (verbose) { std::cout << "ncap  " << buffer.capacity() << ", ncap_expect  " << ncap_expect << std::endl; }
    if (ncap_expect > 0) {
        assert(buffer.size() == nsize_expect);
        assert(buffer.capacity() == ncap_expect);
        if (nsize_expect > 0) {
            assert(!buffer.empty());
            if (nsize_expect == ncap_expect) { assert(buffer.full()); }
            // check raw data
            std::equal(data, data + nsize_expect, buffer.vec().data());
            std::equal(data, data + nsize_expect, buffer.data());
            // check time-ordered indexing
            for (size_t i = 0; i < nsize_expect; ++i) {
                if (verbose) {
                    std::cout << "i " << i << ": buffer[i] " << buffer[i] << ", range[i] " << *(range + i)
                              << ";  buffer.data[i] " << buffer.data()[i] << ", data[i] " << data[i] << std::endl;
                }
                assert(buffer[i] == *(range + i));
            }
            assert(buffer.front() == *(range));
            assert(buffer.back() == *(range + nsize_expect - 1));
        }
        else {
            assert(buffer.empty());
            assert(!buffer.full());
        }
    }
    else {
        assert(buffer.empty());
        assert(buffer.full());
        const size_t size = buffer.size();
        assert(size == 0);
        assert(buffer.capacity() == 0);
    }
}

int main()
{
    using namespace std;
    using namespace nfm;

    const size_t ntest = 12;
    double testvec[ntest];
    std::iota(testvec, testvec + ntest, 0.); // fill with 0..11

    const size_t nbuf = 10;
    PushBackBuffer<double> testbuf(nbuf);

    double wrapvec[nbuf]; // a wrapped version for comparison (here {10, 11, 2, .., 9})
    std::iota(wrapvec, wrapvec + ntest - nbuf, nbuf); // here 10 and 11
    std::iota(wrapvec + ntest - nbuf, wrapvec + nbuf, ntest - nbuf); // 2..9

    const size_t nsmall = 3;

    if (verbose) { cout << "Check the empty buffer:" << endl; }
    assertBuffer(testbuf, nbuf, 0, testvec, testvec);

    // fill in the first nbuf - 1 data elements
    if (verbose) { cout << "Filling " << nbuf - 1 << " elements into ring of cap " << nbuf << ":" << endl; }
    for (size_t i = 0; i < nbuf - 1; ++i) { // leave one element free (i.e. buffer not full)
        testbuf.push_back(testvec[i]);
    }
    assertBuffer(testbuf, nbuf, nbuf - 1, testvec, testvec);

    // increase cap temporarily
    if (verbose) { cout << "Increasing capacity to " << ntest << ":" << endl; }
    testbuf.set_cap(ntest); // calls reserve() ("else" case)
    assertBuffer(testbuf, ntest, nbuf - 1, testvec, testvec);

    // make it one less
    if (verbose) { cout << "Decreasing capacity by 1 again:" << endl; }
    testbuf.set_cap(ntest - 1); // calls "else" case of set_cap()
    assertBuffer(testbuf, ntest - 1, nbuf - 1, testvec, testvec);

    // now decrease it to current size(!) - 1 (dropping oldest element)
    if (verbose) { cout << "Now decreasing capacity to one less than current size(" << nbuf - 1 << "):" << endl; }
    testbuf.set_cap(nbuf - 2); // calls "else if" case of set_cap()
    assertBuffer(testbuf, nbuf - 2, nbuf - 2, testvec+1, testvec+1); // oldest (first) element was dropped

    // now set the buffer back to the original
    if (verbose) { cout << "And back to the original of " << nbuf << ":" << endl; }
    testbuf.push_back(0.5); // push some garbage (should not appear anymore in the assert below)
    for (size_t i = 0; i < nbuf - 2; ++i) { // push through one full round with proper data
        testbuf.push_back(testvec[i]);
    }
    testbuf.set_cap(nbuf); // calls reserve() again, but "if" case
    testbuf.push_back(testvec[nbuf - 2]); // restore the last element
    assertBuffer(testbuf, nbuf, nbuf - 1, testvec, testvec); // we should have the original buffer


    // fill in the rest of the data (which will fill and then cycle the buffer)
    if (verbose) { cout << "Now emplace the remaining " << ntest-nbuf << " elements into the ring:" << endl; }
    for (size_t i = nbuf - 1; i < ntest; ++i) { // now push the rest
        testbuf.emplace_back(testvec[i]); // I don't think emplacing a double means much but we do it for testing
    }
    assertBuffer(testbuf, nbuf, nbuf, testvec + ntest - nbuf, wrapvec);

    // reduce the buffer size significantly
    if (verbose) { cout << "And reduce the cap to " << nsmall << ":" << endl; }
    testbuf.set_cap(nsmall); // keep only most recent data
    assertBuffer(testbuf, nsmall, nsmall, testvec + ntest - nsmall, testvec + ntest - nsmall); // will be ordered again

    if (verbose) { cout << "Try making a copy:" << endl; }
    PushBackBuffer<double> testbuf2(testbuf); // make a copy
    assertBuffer(testbuf2, testbuf.capacity(), testbuf.size(), testbuf.data(), testbuf.data());

    if (verbose) { cout << "Clear the original:" << endl; }
    testbuf.clear();
    assertBuffer(testbuf, nsmall, 0, testvec, testvec);

    if (verbose) { cout << "And check swap:" << endl; }
    swap(testbuf, testbuf2);
    assertBuffer(testbuf, nsmall, nsmall, testvec + ntest - nsmall, testvec + ntest - nsmall); // like above

    return 0;
}
