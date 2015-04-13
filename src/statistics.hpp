#ifndef GI_STATISTICS_HPP
#define GI_STATISTICS_HPP

#include "standard.hpp"
#include "confusion_matrix.hpp"

namespace Spike
{
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

    template <typename T, typename F> void classify(const Standard &r, const T &t, ConfusionMatrix &m, F f)
    {
        const bool mapped = r.l.contains(t.l);

        if (r.id == t.id)
        {
            if (mapped)
            {
                if (f(t))
                {
                    m.tp++;
                }
                else
                {
                    m.fp++;
                }
            }
            else
            {
                m.fp++;
            }
        }
        else
        {
            if (mapped)
            {
                m.fn++;
            }
            else
            {
                m.tn++;
            }
        }
    }
}

#endif