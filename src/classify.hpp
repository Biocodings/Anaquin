#ifndef GI_CLASSIFY_HPP
#define GI_CLASSIFY_HPP

#include "standard.hpp"
#include <ss/ml/classify.hpp>

namespace Spike
{
    struct Confusion : public SS::Confusion
    {
        inline Percentage sp() const
        {
            return (tp + fp) ? tp / (tp + fp) : NAN;
        }
    };
    
    inline bool tfp(bool cond, Confusion *m1, Confusion *m2 = NULL)
    {
        if (cond)
        {
            if (m1) { m1->tp++; }
            if (m2) { m2->tp++; }
        }
        else
        {
            if (m1) { m1->fp++; }
            if (m2) { m2->fp++; }
        }

        return cond;
    }

    template <typename Iter, typename T> bool find(const Iter &iter, const T &t)
    {
        for (auto i: iter)
        {
            if (i.l.contains(t.l))
            {
                return true;
            }
        }

        return false;
    }

    /*
     * Specalised binary-classification for the project. Negativity isn't required because in a typical
     * experiment the dilution would be so low that true-negative dominates false-positive.
     */
    
    template <typename T, typename Stats, typename Positive>
    void classify(Stats &stats, const T &t, Positive p)
    {
        static const Standard &r = Standard::instance();

        SS::classify(stats.m_base, t,
            [&](const T &)  // Classifier
            {
                return (t.id == r.id && r.l.contains(t.l));
            },
            [&](const T &t) // Positive (reference detected)
            {
                stats.n++;
                stats.nr++;
                return p(t);
            },
            [&](const T &t) // Negative (query detected)
            {
                stats.n++;
                stats.nq++;
                return false;
            });
       }
}

#endif