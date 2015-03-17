#ifndef AS_STATISTICS_HPP
#define AS_STATISTICS_HPP

#include "ConfusionMatrix.hpp"

template <typename T, typename Iter> bool contains(const Iter &iter, const T &t)
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
 * Conduce a binary classification test. Refer to http://en.wikipedia.org/wiki/Sensitivity_and_specificity for more details.
 */

template <typename Iter, typename R, typename T> void binaryClassify(const Iter &iter, const R &r, const T &t, ConfusionMatrix &m)
{
    assert(!iter.empty());
    
    if (t.chromo == r.id)
    {
        if (contains(iter, t))
        {
            if (r.l.contains(t.l))
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
        if (contains(iter, t))
        {
            m.fn++;
        }
        else
        {
            m.tn++;
        }
    }
}

#endif