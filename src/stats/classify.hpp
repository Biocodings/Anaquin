#ifndef CLASSIFY_HPP
#define CLASSIFY_HPP

#include "stats/limit.hpp"
#include "data/standard.hpp"

namespace Anaquin
{
    class Confusion
    {
        public:
        
            inline Counts &tp() const { return _tp; }
            inline Counts &fp() const { return _fp; }
            inline Counts &tn() const { return _tn; }
            inline Counts &fn() const { return _fn; }

            inline Counts &nr() const { return _nr; }
            inline Counts &nq() const { return _nq; }

            inline Counts n() const
            {
                return _fp + _tp + _fn + _tn;
            }
        
            // Sensitivity, metrics for positive classification
            inline Proportion sn() const
            {
                assert(_nr && _nr >= _tp);

                // Adjust for fn... Refer to wikipedia for details
                _fn = _nr - _tp;

                return (_tp + _fn) ? static_cast<Proportion>(_tp) / (_tp + _fn) : NAN;
            }

            // Specificity, metrics for negative classification
            inline Proportion sp() const
            {
                assert(_nr && _nr >= _tp);
                return (_tn + _fp) ? static_cast<Proportion>(_tn) / (_tn + _fp) : NAN;
            }

            // Precision, metrics for accuracy
            inline Proportion pc() const
            {
                assert(_nr && _nr >= _tp);
                return ((tp() + fp()) && fp() != n()) ? static_cast<Proportion>(tp()) / (tp() + fp()) : NAN;
            }

        private:
        
            mutable Counts _fp = 0;
            mutable Counts _tp = 0;
            mutable Counts _fn = 0;
            mutable Counts _tn = 0;
            mutable Counts _nq = 0;
            mutable Counts _nr = 0;
    };

    /*
     * Overall performance for a metric
     */
    
    struct Performance
    {
        Confusion m;

        // Absolute detection limit
        Limit limit;

        // Histogram of distribution
        SequinHist hist;
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
                m.nq()++;
                m.fp()++;
            }
            else
            {
                m.nq()++;
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