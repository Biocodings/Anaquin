#ifndef SS_F_DIST_HPP
#define SS_F_DIST_HPP

#include <ss/prob.hpp>
#include <ss/stats.hpp>
#include <boost/math/distributions/fisher_f.hpp>

namespace SS
{
    template <typename T> struct FDist
    {
        FDist(DF n, DF d);
    
        T mean() const;
        T variance() const;
        P cdf(T x) const;
        
        DF n, d;
    };
    
    template <typename T> FDist<T>::FDist(DF n, DF d) : n(n), d(d)
    {
        // Empty Implementation
    }
    
    template <typename T> P FDist<T>::cdf(T x) const
    {
        boost::math::fisher_f_distribution<T> c(n, d);
        return P(boost::math::cdf(c, x));
    }

    template <typename T> T FDist<T>::mean() const
    {
        boost::math::fisher_f_distribution<T> f(n, d);
        return boost::math::mean(f);
    }
}

#endif
