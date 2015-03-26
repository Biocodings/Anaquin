#ifndef GI_STATISTICS_HPP
#define GI_STATISTICS_HPP

#include "confusion_matrix.hpp"

template <typename T, typename Iter> bool contains_(const Iter &iter, const T &t)
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

template <typename Iter, typename R, typename T> void classify(const Iter &iter, const R &r, const T &t, ConfusionMatrix &m)
{
    assert(!iter.empty());
    
    // Positive if the query is classified with the references
    if (r.id == t.id)
    {
        if (contains_(iter, t))
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
        if (contains_(iter, t))
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