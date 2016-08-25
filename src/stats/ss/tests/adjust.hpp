#ifndef SS_ADJUST_HPP
#define SS_ADJUST_HPP

#include <ss/internal/rank.hpp>

namespace SS
{
    /*
     * Implements adjustment for multiple testing
     *
     *    - FDR
     */

    template <typename T> T adjustFDR(const T &p)
    {
        std::vector<Index> o;

        // Ordered index for p
        o.resize(p.size());
        
        Internal::order(p, o,  false);
        const auto ro = Internal::getPerm(o);

        auto sorted = p;
        std::sort(sorted.rbegin(), sorted.rend());

        const auto n = p.size();
        auto i = n;

        for (auto &j : sorted)
        {
            j *= static_cast<Real>(n)/i--;
        }
        
        auto runs = sorted;
        Internal::cummin(sorted, runs);
        
        for (auto &j : runs)
        {
            j = std::min(j, 1.0);
        }

        return Internal::applyPerm(runs, ro);
    }
}

#endif