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
                /*
                 * It's a better sensitivity measure if it has smaller counts while still being
                 * detectable. Choose the smaller abundance if there's a tie.
                 */

                if (iter->second)
                {
                    if (iter->second < s.counts || (iter->second == s.counts && m.at(iter->first).abund() < s.abund))
                    {
                        s.id     = iter->first;
                        s.counts = iter->second;
                        //s.abund  = m.at(s.id).abund(); // TODO: Reactive this
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
