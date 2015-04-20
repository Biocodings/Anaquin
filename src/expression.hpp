#ifndef GI_EXPRESSION_HPP
#define GI_EXPRESSION_HPP

#include <map>
#include <iostream>
#include "standard.hpp"
#include "sensitivity.hpp"

namespace Spike
{
    struct Expression
    {
        template <typename T> struct ExpressionResults
        {
            T limit_key;

            // The value of the limit of sensitivity
            unsigned limit_count;

            inline Sensitivity sens() const
            {
                static const auto &r = Standard::instance();

                Sensitivity s;
                
                s.id     = limit_key;
                s.counts = limit_count;
                s.exp    = limit_count ? r.seqs_iA.at(limit_key).raw +
                                         r.seqs_iA.at(limit_key).raw: NAN;
                return s;
            }
        };

        template <typename T> static void print(const std::map<T, unsigned> &m)
        {
            for (auto iter = m.begin(); iter != m.end(); iter++)
            {
                std::cout << iter->first << "  " << iter->second << std::endl;
            }
        }

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