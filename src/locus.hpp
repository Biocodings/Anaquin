#ifndef GI_LOCUS_HPP
#define GI_LOCUS_HPP

#include <assert.h>
#include "types.hpp"

namespace Spike
{
    struct Locus
    {
        Locus(BasePair start = 0, BasePair end = 0) : start(start), end(end) {}
        
        inline void set(BasePair start, BasePair end)
        {
            this->start = start; this->end = end;
            assert(this->end >= this->start);
        }

        inline BasePair length() const { return (end - start + 1); }

        // Returns the base-pairs that are overlapped with the given reference
        inline BasePair overlap(const Locus &l) const
        {
            if (l.start > end || start > l.end)
            {
                return 0;
            }
            else if (start <= l.start && end >= l.end)
            {
                return l.length();
            }
            else if (start >= l.start && end <= l.end)
            {
                return length();
            }
            else if (end >= l.end)
            {
                return l.end - start + 1;
            }
            else
            {
                return end - l.start + 1;
            }
        }

        inline bool contains(const Locus &q) const
        {
            return (q.start >= start && q.end <= end);
        }

        inline Locus operator+(const Locus &l) const
        {
            return Locus(std::min(start, l.start), std::max(end, l.end));
        }

        bool operator==(const Locus &x) const
        {
            return start == x.start && end == x.end;
        }

        BasePair start, end;
    };    
}

#endif