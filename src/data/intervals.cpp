#include "data/intervals.hpp"

using namespace Anaquin;

Interval::Stats Intervals::stats() const
{
    Interval::Stats stats;
    
    for (const auto &i : _inters)
    {
        const auto s = i.second.stats();
        
        stats.sums     += s.sums;
        stats.length   += s.length;
        stats.nonZeros += s.nonZeros;
        stats.zeros    += s.zeros;
        stats.min       = std::min(stats.min, s.min);
        stats.max       = std::max(stats.max, s.max);
        
        for (const auto &j : s.hist)
        {
            stats.hist[j.first] += j.second;
        }
    }

    stats.mean = stats.sums / stats.length;

    return stats;
}