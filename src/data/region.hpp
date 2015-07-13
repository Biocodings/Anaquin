#ifndef GI_REGION_HPP
#define GI_REGION_HPP

#include <vector>
#include <numeric>
#include "data/locus.hpp"

namespace Anaquin
{
    struct Region : public std::vector<Locus>
    {
        inline Locus region() const
        {
            auto end   = std::numeric_limits<BasePair>::min();
            auto start = std::numeric_limits<BasePair>::max();

            for (const auto &l : *this)
            {
                end   = std::max(end, l.end);
                start = std::min(start, l.start);
            }

            return Locus(start, end);
        }
    };
}

#endif
