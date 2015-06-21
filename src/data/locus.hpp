#ifndef GI_LOCUS_HPP
#define GI_LOCUS_HPP

#include <set>
#include <map>
#include <list>
#include <vector>
#include <assert.h>
#include <algorithm>
#include "data/types.hpp"

namespace Spike
{
    struct Locus
    {
        Locus(const Locus &l1, const Locus &l2)
        {
            end   = std::max(l1.end,   l2.end);
            start = std::min(l1.start, l2.start);
        }

        Locus(BasePair start = 0, BasePair end = 0) : start(start), end(end) {}

        template <typename Iter, typename F> static Locus expand(const Iter &iter, F f)
        {
            Locus l;
            
            l.end   = std::numeric_limits<BasePair>::min();
            l.start = std::numeric_limits<BasePair>::max();

            for (const auto &i : iter)
            {
                if (f(i))
                {
                    l.end   = std::max(l.end, i.l.end);
                    l.start = std::min(l.start, i.l.start);
                }
            }

            return l;
        };
        
        template <typename T> static bool overlap(const std::vector<T> &ls)
        {
            for (auto i = 0; i < ls.size(); i++)
            {
                for (auto j = i + 1; j < ls.size(); j++)
                {
                    if (ls[i].overlap(ls[j]))
                    {
                        return true;
                    }
                }
            }

            return false;
        }

        /*
         * Merge a list of objects according to the locus. The type of the object must be
         * castable to Locus.
         */
        
        template <typename T, typename R, template <typename, typename = std::allocator<T>> class Inputs>
                static std::vector<R> merge(const Inputs<T> &x)
        {
            std::vector<R> sorted;

            // std::copy doesn't allow implicit type conversion...
            for (const auto &i : x)
            {
                sorted.push_back(i);
            }

            std::sort(sorted.begin(), sorted.end(), [&](const Locus &l1, const Locus &l2)
            {
                return (l1.start < l2.start) || (l1.start == l2.start && l1.end < l2.end);
            });
            
            std::vector<R> merged;
            
            for (auto i = 0; i < sorted.size();)
            {
                // We'll need the index for skipping i
                auto j = i + 1;

                auto super = sorted[i];
                
                // Look forward until the end of an overlap region
                for (; j < sorted.size(); j++)
                {
                    if (sorted[j].overlap(super))
                    {
                        super += sorted[j];
                    }
                    else
                    {
                        break;
                    }
                }

                // Construct the super-loci for the region
                merged.push_back(super);

                i = j;
            }

            return merged;
        }
        
        inline BasePair length() const { return (end - start + 1); }

        inline BasePair overlap(const Locus &l) const
        {
            if (l.start > end || start > l.end)
            {
                return 0;
            }
            else if (start <= l.start && end >= l.end)
            {
                return l.length();
            }
            else if (start >= l.start && end <= l.end)
            {
                return length();
            }
            else if (end >= l.end)
            {
                return l.end - start + 1;
            }
            else
            {
                return end - l.start + 1;
            }
        }

        inline bool contains(const Locus &q) const
        {
            return (q.start >= start && q.end <= end);
        }

        inline void operator+=(const Locus &l)
        {
            this->start = std::min(start, l.start);
            this->end   = std::max(end, l.end);
        }

        inline bool operator!=(const Locus &l) const { return !operator==(l); }
        inline bool operator==(const Locus &l) const { return start == l.start && end == l.end; }

        inline Locus operator+(const Locus &l)  const
        {
            return Locus(std::min(start, l.start), std::max(end, l.end));
        }

        inline bool operator<(const Locus &l)  const
        {
            return start < l.start || (start == l.start && end < l.end);
        }

        BasePair start, end;
    };

    struct RNALocus : public Locus
    {
        RNALocus(const std::string &gID, const Locus &l) : gID(gID), Locus(l) {}

        std::string gID;
    };
}

#endif
