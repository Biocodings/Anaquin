#ifndef CLASSIFY_HPP
#define CLASSIFY_HPP

#include "stats/limit.hpp"
#include "data/standard.hpp"
#include <ss/data/confusion.hpp>

namespace Anaquin
{
    struct Confusion : public SS::Confusion
    {
        inline Counts &tn() const { throw std::runtime_error("TN is unsupported"); }

        // Sensitivity, metrics for positive classification
        inline Percentage sn() const
        {
            assert(nr && nr >= _tp);

            // Adjust for fn... Refer to the wikipedia for more details
            _fn = nr - _tp;

            return SS::Confusion::sn();
        }
        
        /*
         * The usual formula: tn / (tn + fp) would not work because we don't know tn.
         * Furthermore fp would have been dominated by tn. The formula below is
         * consistent to cufflink's recommendation. Technically, we're not calculating
         * specificity but prediction accuracy.
         *
         * Reference: https://www.biostars.org/p/138438/#148051
         */

        inline Percentage sp() const
        {
            return ((tp() + fp()) && fp() != n()) ? static_cast<Percentage>(tp()) / (tp() + fp()) : NAN;
        }
        
        inline Percentage accuracy() const
        {
            return ((tp() + fp()) && fp() != n()) ? static_cast<Percentage>(tp()) / (tp() + fp()) : NAN;
        }

        // Counts of queries
        Counts nq = 0;
        
        // Counts of references
        Counts nr = 0;
    };

    /*
     * Overall performance for a metric
     */
    
    struct Performance
    {
        Confusion m;

        // Hard detection limit
        Limit hl;

        // Soft detection limit
        Limit sl;

        // Histogram of distribution
        Hist h;

        // Calculate the references from the distribution
        inline void inferRefFromHist()
        {
            for (const auto &i : h)
            {
                if (i.second == 0)
                {
                    m.nr++;
                }
                else
                {
                    m.nr += i.second;
                }
            }
        }
    };

    /*
     * Given references, this function computes total overalap of a query to the references.
     * For example:
     *
     *    r = (1,10),(12,20)  q = (15,21)
     *
     * The overlap consists of (15,20).
     */

    template <typename Iter, typename T, typename C> Base countOverlaps(const Iter &rs, const T &q, C &c)
    {
        Base n = 0;

        for (const auto &x : rs)
        {
            if (static_cast<Locus>(q).overlap(x))
            {
                c.at(x.gID)++;
                n += static_cast<Locus>(q).overlap(x);
            }
        }

        return n;
    }

    template <typename Iter, typename T> const typename Iter::value_type * find
        (const Iter &iter, const T &t, MatchRule rule)
    {
        for (const auto &i : iter)
        {
            const bool matched = (rule == Exact    && i.l == static_cast<Locus>(t)) ||
                                 (rule == Contains && i.l.contains(static_cast<Locus>(t)));
            if (matched)
            {
                return &i;
            }
        }

        return NULL;
    }

    // Similar to find(), works on a std::map
    template <typename Key, typename Value, typename T> const typename std::map<Key, Value>::mapped_type * find
        (const std::map<Key, Value> &m, const T &t, MatchRule rule)
    {
        for (const auto &i: m)
        {
            const auto matched = (rule == Exact    && i.second.l == t.l) ||
                                 (rule == Contains && i.second.l.contains(t.l));
            if (matched)
            {
                return &i.second;
            }
        }

        return NULL;
    }

    enum ClassifyResult
    {
        Ignore   = -1,
        Negative = 0,
        Positive = 1,
    };
    
    template <typename T, typename Classifer> bool classify(Confusion &m, const T &t, Classifer c)
    {
        const auto r = c(t);
        
        if ((void *) r != (void *) ClassifyResult::Ignore)
        {
            if ((void *) r == (void *) ClassifyResult::Negative)
            {
                m.nq++;
                m.fp()++;
            }
            else
            {
                m.nq++;
                m.tp()++;
                return true;
            }
        }

        return false;
    }
    
    template <typename T> bool classifyTP(Confusion &m, const T &t)
    {
        return classify(m, t, [&](const T &t)
        {
            return true;
        });
    }

    template <typename T> bool classifyFP(Confusion &m, const T &t)
    {
        return classify(m, t, [&](const T &t)
        {
            return false;
        });
    }
}

#endif