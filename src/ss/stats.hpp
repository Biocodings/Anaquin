#ifndef SS_STATS_HPP
#define SS_STATS_HPP

#include <map>
#include <math.h>
#include <numeric>
#include <limits.h>
#include <stdexcept>
#include <ss/types.hpp>

namespace SS
{
    namespace Internal
    {
        template <class Iter> Correlation kendall(Iter x, Iter y)
        {
            const auto nx = count(x);
            const auto ny = count(y);
            
            if (nx != ny)
            {
                throw std::runtime_error("Incompatible dimensions");
            }
            
            Counts n = 0;
            
            for (auto i = 0; i < nx; i++)
            {
                for (auto j = 0; j < i; j++)
                {
                    n = n + sign(x[i] - x[j]) * sign(y[i] - y[j]);
                }
            }
            
            return (n / (0.5 * nx * (nx - 1)));
        }
        
        template <typename Iter> typename Iter::value_type sum(const Iter &x)
        {
            return std::accumulate(x.begin(), x.end(), 0.0);
        }
        
        template <typename Iter> Counts count(const Iter &x)
        {
            return static_cast<Counts>(std::distance(x.begin(), x.end()));
        }
        
        template <typename T> int sign(T val)
        {
            return (T(0) < val) - (val < T(0));
        }

        template <class Iter> typename Iter::value_type var(const Iter &x)
        {
            typedef typename Iter::value_type T;
            
            const auto sum = std::accumulate(x.begin(), x.end(), 0.0);
            const auto mean = sum / x.size();
            
            std::vector<T> diff(x.size());
            std::transform(x.begin(), x.end(), diff.begin(),
                           std::bind2nd(std::minus<T>(), mean));

            const auto sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
            return sq_sum / (x.size() - 1);
        }

        template <typename Iter> typename Iter::value_type sd__(Iter x)
        {
            return sqrt(var(x));
        }

        template <class Iter> typename Iter::value_type cov(const Iter &x, const Iter &y)
        {
            const auto sum_x = std::accumulate(std::begin(x), std::end(x), 0.0);
            const auto sum_y = std::accumulate(std::begin(y), std::end(y), 0.0);

            const auto mx =  sum_x / x.size();
            const auto my =  sum_y / y.size();

            auto accum = 0.0;

            for (auto i = 0; i < x.size(); i++)
            {
                accum += (x.at(i) - mx) * (y.at(i) - my);
            }
            
            return accum / (x.size() - 1);
        }

        template <typename Iter> Correlation pearson(const Iter &x, const Iter &y)
        {
            return (cov(x, y) / (sd__(x) * sd__(y)));
        }
    }

    /*
     * Generic function for the (trimmed) arithmetic mean.
     */

    template <typename Iter> Mean mean(const Iter &x)
    {
        return Internal::sum(x) / Internal::count(x);
    }

    /*
     * This function computes the standard deviation of the values in x
     */
    
    template <typename Iter> SD sd(const Iter &x)
    {
        return Internal::sd__(x);
    }
    
    /*
     * This function computes the variance of the values in x
     */
    
    template <typename Iter> SD var(const Iter &x)
    {
        return Internal::var(x);
    }

    /*
     * This function computes the correlation between x and y.
     */
    
    template <typename Iter> Correlation cor(const Iter &x, const Iter &y)
    {
        return Internal::pearson(x, y);
    }

    /*
     * This function computes the covariance between x and y.
     */
    
    template <typename Iter> typename Iter::value_type cov(const Iter &x, const Iter &y)
    {
        return Internal::cov(x, y);
    }
}

#endif