#ifndef SF_CHI_SQUARE_HPP
#define SF_CHI_SQUARE_HPP

#include <sstats/stats.hpp>
#include <sstats/types.hpp>
#include <sstats/probability.hpp>
#include <boost/math/distributions/chi_squared.hpp>

namespace SS
{
    namespace Internal
    {
        template <typename T> struct ChiSquare
        {
            ChiSquare(DegreeFreedom df) : df(df) {}
            
            Probability<T> cdf(double x) const;
            
            DegreeFreedom df;
        };

        template <typename T> Probability<T> ChiSquare<T>::cdf(double x) const
        {
            boost::math::chi_squared_distribution<T> c(df);
            return Probability<T>(boost::math::cdf(c, x));
        }
    }

    typedef Internal::ChiSquare<Real> ChiSquare;
}

#endif