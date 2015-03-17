#ifndef AS_LOCUS_HPP
#define AS_LOCUS_HPP

#include <assert.h>
#include "Types.hpp"

struct Locus
{
    Locus(BasePair start = 0, BasePair end = 0) : start(start), end(end) {}
    ~Locus() {}

    inline void update(BasePair start, BasePair end)
    {
        this->start = start; this->end = end;
        assert(this->end > start);
    }

    inline BasePair length() const { return (end - start); }

    inline bool contains(const Locus &q) const
    {
        return (q.start >= start && q.end <= end);
    }

    BasePair start, end;
};

#endif