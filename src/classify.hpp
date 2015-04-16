#ifndef GI_CLASSIFY_HPP
#define GI_CLASSIFY_HPP

#include "standard.hpp"
#include "confusion.hpp"

namespace Spike
{
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

    inline bool tfn(bool cond, Confusion *m1, Confusion *m2 = NULL)
    {
        if (cond)
        {
            if (m1) { m1->tn++; }
            if (m2) { m2->tn++; }
        }
        else
        {
            if (m1) { m1->fn++; }
            if (m2) { m2->fn++; }
        }

        return cond;
    }

    template <typename Iter, typename F> void extractIntrons(const Iter &exons, F f)
    {
        Feature in;
        
        for (auto i = 0; i < exons.size(); i++)
        {
            if (i)
            {
                in.l = Locus(exons[i - 1].l.end, exons[i].l.start);                
                f(exons[i - 1], exons[i], in);
            }
        }
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
     * Define a generic algorithm for experimental classification. The caller is
     * expected to provide a functor for positive and negative match.
     */
    
    template <typename T, typename Stats, typename Positive, typename Negative>
    void classify(const Standard &r, Stats &stats, const T &t, Positive p, Negative n)
    {
        if (t.id == r.id)
        {
            stats.nr++;
        }
        else
        {
            stats.nq++;
        }

        stats.n++;

        // Whether it's been mapped to the reference
        const bool mapped = r.l.contains(t.l);

        if (r.id == t.id && mapped)
        {
            p();
        }
        else
        {
            n(mapped);
        }
    }
}

#endif