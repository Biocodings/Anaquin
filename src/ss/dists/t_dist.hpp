#ifndef SS_T_DIST_HPP
#define SS_T_DIST_HPP

#include <ss/prob.hpp>
#include <ss/types.hpp>
#include <boost/math/distributions/students_t.hpp>

namespace SS
{
    template <typename T> struct TDist
    {
        TDist(DF df);
        
        P cdf(double x) const;
        
        // The distribution is specified by it's degree-of-freedom
        DF df;
    };
    
    template <typename T> TDist<T>::TDist(DF df) : df(df) {}
    
    template <typename T> P TDist<T>::cdf(double x) const
    {
        boost::math::students_t_distribution<T> t(df);
        return boost::math::cdf(t, x);
    }
    
    namespace R
    {
        /*
         * Density, distribution function, quantile function and random generation
         * for the t distribution with df degrees of freedom (and optional non-centrality
         * parameter ncp).
         */

        inline P pt(double x, DF df)
        {
            return TDist<double>(df).cdf(x);
        }
    }
}

#endif