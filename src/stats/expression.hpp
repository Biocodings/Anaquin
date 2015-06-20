#ifndef GI_EXPRESSION_HPP
#define GI_EXPRESSION_HPP

#include <map>
#include <limits>
#include <math.h>
#include <iostream>
#include <assert.h>
#include "sensitivity.hpp"

namespace Spike
{
    struct Expression
    {
        template <typename Map> static void print(const Map &m)
        {
            for (auto iter = m.begin(); iter != m.end(); iter++)
            {
                std::cout << iter->first << "  " << iter->second << std::endl;
            }
        }

        template <typename T, typename ID, typename S> static Sensitivity
                analyze(const std::map<T, Counts> &c, const std::map<ID, S> &m)
        {
            Sensitivity s;

            // The lowest count must be zero because it can't be negative
            s.counts = std::numeric_limits<unsigned>::max();

            for (auto iter = c.begin(); iter != c.end(); iter++)
            {
                const auto counts = iter->second;

                /*
                 * Is this sequin detectable? If it's detectable, what about the concentration?
                 * By definition, the detection limit is defined as the smallest abundant while
                 * still being detected.
                 */

                if (counts)
                {
                    if (counts < s.counts || (counts == s.counts && m.at(iter->first).abund() < s.abund))
                    {
                        s.id     = iter->first;
                        s.counts = counts;
                        s.abund  = m.at(s.id).abund();
                    }
                }
            }

            if (s.counts == std::numeric_limits<unsigned>::max())
            {
                s.counts = 0;
            }

            return s;
        }
    };
}

#endif
