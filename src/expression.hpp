#ifndef GI_EXPRESSION_HPP
#define GI_EXPRESSION_HPP

#include <map>

namespace Spike
{
    struct Expression
    {
        template <typename T> struct ExpressionResults
        {
            T limit_key;

            // The value of the limit of sensitivity
            unsigned limit_count;
        };

        template <typename T> static ExpressionResults<T> analyze(const std::map<T, unsigned> &t)
        {
            ExpressionResults<T> r;

            if (!t.size()) { return r; }

            // The lowest count must be zero because it can't be negative
            r.limit_count = std::numeric_limits<unsigned>::max();

            for (auto iter = t.begin(); iter != t.end(); iter++)
            {
                if (iter->second && iter->second < r.limit_count)
                {
                    r.limit_key = iter->first;
                    r.limit_count = iter->second;
                }
            }

            if (r.limit_count == std::numeric_limits<unsigned>::max())
            {
                // There's nothing that has a single count
                r.limit_count = 0;
            }
            
            return r;
        }
    };
}

#endif