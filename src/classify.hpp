#ifndef GI_CLASSIFY_HPP
#define GI_CLASSIFY_HPP

#include "standard.hpp"
#include "confusion_matrix.hpp"

namespace Spike
{
    void verifytPositive(bool cond, ConfusionMatrix *m1, ConfusionMatrix *m2 = NULL)
    {
        if (cond)
        {
            if (m1) { m1->tp++; }
            if (m2) { m2->tp++; }
        }
        else
        {
            if (m1) { m1->tn++; }
            if (m2) { m2->tn++; }
        }
    }

    void verifyNegative(bool cond, ConfusionMatrix *m1, ConfusionMatrix *m2 = NULL)
    {
        if (cond)
        {
            if (m1) { m1->fp++; }
            if (m2) { m2->fp++; }
        }
        else
        {
            if (m1) { m1->fn++; }
            if (m2) { m2->fn++; }
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

    template <typename T, typename Positive, typename Negative>
    void classify(const Standard &r, const T &t, Positive p, Negative n)
    {
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