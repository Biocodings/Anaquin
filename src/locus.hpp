#ifndef GI_LOCUS_HPP
#define GI_LOCUS_HPP

#include <assert.h>
#include "types.hpp"

namespace Spike
{
    struct Locus
    {
        Locus(BasePair start = 0, BasePair end = 0) : start(start), end(end) {}
        ~Locus() {}
        
        inline void set(BasePair start, BasePair end)
        {
            this->start = start; this->end = end;
            assert(this->end > this->start);
        }

        inline BasePair length() const { return (end - start + 1); }
        
        inline bool contains(const Locus &q) const
        {
            return (q.start >= start && q.end <= end);
        }
        
        BasePair start, end;
    };    
}

#endif