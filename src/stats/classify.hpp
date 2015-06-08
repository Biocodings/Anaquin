#ifndef GI_CLASSIFY_HPP
#define GI_CLASSIFY_HPP

#include "sensitivity.hpp"
#include "data/standard.hpp"
#include <ss/data/confusion.hpp>

namespace Spike
{
    struct Confusion : public SS::Confusion
    {
        inline Counts &tn() const { throw std::runtime_error("tn() is unsupported"); }
        
        inline Percentage sn() const
        {
            assert(nr && nr >= _tp);
            
            // Adjust for fn... Refer to the wikipedia for more details
            _fn = nr - _tp;
            
            return SS::Confusion::sn();
        }
        
        /*
         * The usual formula: tn / (tn + fp) would not work. We don't know
         * tn, furthermore fp would have been dominated by tn. The formula
         * below is consistent to cufflink's recommendation. Technically,
         * we're not calculating specificity but positive predication value.
         */
        
        inline Percentage sp() const
        {
            return ((tp() + fp()) && fp() != n()) ? static_cast<Percentage>(tp()) / (tp() + fp()) : NAN;
        }

        Counts nq = 0;
        Counts nr = 0;
    };

    struct Performance
    {
        Confusion m;
        Sensitivity s;
    };

    inline bool tfp(bool cond, Confusion *m1, Confusion *m2 = NULL)
    {
        if (cond)
        {
            if (m1) { m1->tp()++; }
            if (m2) { m2->tp()++; }
        }
        else
        {
            if (m1) { m1->fp()++; }
            if (m2) { m2->fp()++; }
        }

        return cond;
    }

    template <typename Iter, typename T> BasePair countOverlaps(const Iter &r, const T &t, std::map<std::string, Counts> &c)
    {
        BasePair n = 0;

        for (const auto &i : r)
        {
            if (static_cast<Locus>(t).overlap(i))
            {
                c.at(i.gID)++;
                n += static_cast<Locus>(t).overlap(i);
            }
        }
        
        return n;
    }

    template <typename Iter, typename T> BasePair countOverlaps_map(const Iter &map, const T &t)
    {
        BasePair n = 0;
        
        for (auto &i : map)
        {
            n += i.second.l.overlap(t.l);
        }

        return n;
    }

    enum MatchRule
    {
        Exact,
        Contains,
    };

    template <typename Iter, typename T> const typename Iter::value_type * find(const Iter &iter, const T &t, MatchRule rule)
    {
        for (const auto &i : iter)
        {
            const bool matched = (rule == Exact    && i.l == t.l) ||
                                 (rule == Contains && i.l.contains(t.l));
            if (matched)
            {
                return &i;
            }
        }

        return NULL;
    }

    template <typename Iter, typename T> bool find_map(const Iter &map, const T &t, MatchRule rule)
    {
        for (auto i: map)
        {
            const bool matched = (rule == Exact    && i.second.l == t.l) ||
                                 (rule == Contains && i.second.l.contains(t.l));
            if (matched)
            {
                return true;
            }
        }

        return false;
    }

    enum ClassifyResult
    {
        Ignore = -1,
        Negative = 0,
        Positive = 1,
    };
    
    template <typename T, typename Classifer> bool classify(Confusion &m, const T &t, Classifer c)
    {
        const auto &s = Standard::instance();

        if (t.id == s.id && s.l.contains(t.l))
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
        }

        return false;
    }
}

#endif