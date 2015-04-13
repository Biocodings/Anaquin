#ifndef GI_ABUNDANCE_HPP
#define GI_ABUNDANCE_HPP

#include <map>

namespace Spike
{
    struct Abdunance
    {
        template <typename T1, typename T2> struct AbdunanceResults
        {
            T1 min_k;
            T2 min_v;
        };
        
        template <typename T1, typename T2> static AbdunanceResults<T1, T2> analyze(const std::map<T1, T2> &t)
        {
            AbdunanceResults<T1, T2> r;

            if (!t.size()) { return r; }
            
            r.min_v = std::numeric_limits<T2>::max();
            
            for (auto iter = t.begin(); iter != t.end(); iter++)
            {
                if ( iter->second && iter->second < r.min_v)
                {
                    r.min_k = iter->first;
                    r.min_v = iter->second;
                }
            }

            assert(r.min_v);
            return r;
        }
    };
}

#endif