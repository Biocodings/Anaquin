#ifndef GI_EXPRESSION_HPP
#define GI_EXPRESSION_HPP

#include <math.h>
#include <iostream>
#include "standard.hpp"
#include "sensitivity.hpp"

namespace Spike
{
    struct Expression
    {
        template <typename T> struct ExpressionResults
        {
            // Name of a sequin or a gene
            T key;

            // The counts for the sensitivity
            Counts counts;

            template <typename ID, typename S> Sensitivity sens(const std::map<ID, S> &m) const
            {
                Sensitivity s;
                
                s.id     = key;
                s.counts = counts;
                s.abund  = counts ? m.at(key).abund(false) : NAN;

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

        template <typename T> static ExpressionResults<T> analyze(const std::map<T, Counts> &c)
        {
            ExpressionResults<T> r;

            if (!c.size()) { return r; }

            // The lowest count must be zero because it can't be negative
            r.counts = std::numeric_limits<unsigned>::max();

            for (auto iter = c.begin(); iter != c.end(); iter++)
            {
                if (iter->second && iter->second < r.counts)
                {
                    r.key = iter->first;
                    r.counts = iter->second;
                }
            }

            if (r.counts == std::numeric_limits<unsigned>::max())
            {
                r.counts = 0;
            }
            
            return r;
        }
    };
}

#endif