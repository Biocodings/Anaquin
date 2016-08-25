#ifndef SS_INTERNAL_STATS_HPP
#define SS_INTERNAL_STATS_HPP

#include <numeric>
#include <algorithm>
#include <ss/data/data.hpp>
#include <ss/internal/rank.hpp>

namespace SS
{
    namespace Internal
    {
        template <typename T> typename T::value_type quantile(const T &x, P probs)
        {
            const auto n  = x.size();
            const auto id = (n - 1) * probs;
            const auto lo = floor(id);
            const auto hi = ceil(id);
            const auto qs = x[lo];
            const auto h  = (id - lo);
            
            return (1.0 - h) * qs + h * x[hi];
        }
        
        template <typename T> typename T::value_type sum(const T &x)
        {
            return std::accumulate(x.begin(), x.end(), 0.0);
        }
        
        template <typename T> Counts count(const T &x)
        {
            return static_cast<Counts>(std::distance(x.begin(), x.end()));
        }
        
        template <typename T> typename T::value_type cov(const T &x, const T &y)
        {
            const auto sum_x = std::accumulate(std::begin(x), std::end(x), 0.0);
            const auto sum_y = std::accumulate(std::begin(y), std::end(y), 0.0);
            
            const auto mx =  sum_x / x.size();
            const auto my =  sum_y / y.size();
            
            auto accum = 0.0;
            
            auto i = x.begin();
            auto j = y.begin();
            
            for (; i != x.end(); i++, j++)
            {
                accum += (*i - mx) * (*j - my);
            }
            
            return accum / (x.size() - 1.0);
        }

        template <class T> typename T::value_type getVariance(const T &x)
        {
            typedef typename T::value_type Value;
            
            const auto sum  = std::accumulate(x.begin(), x.end(), 0.0);
            const auto mean = sum / x.size();
            
            std::vector<Value> diff(x.size());
            std::transform(x.begin(), x.end(), diff.begin(), std::bind2nd(std::minus<Value>(), mean));
            
            const auto sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
            return sq_sum / (x.size() - 1);
        }
        
        // Computes the standard deviation
        template <typename T> typename T::value_type getSD(const T &x)
        {
            return sqrt(getVariance(x));
        }

        template <typename T> Real getCovariance(const T &x, const T &y)
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
        
        template <typename T> Real corrPearson(const T &x, const T &y)
        {
            return getCovariance(x, y) / (getSD(x) * getSD(y));
        }
        
        template <typename T> Real corrSpearman(const T &x, const T &y)
        {
            std::vector<Rank> rx;
            std::vector<Rank> ry;
            
            rx.resize(x.size());
            ry.resize(y.size());
            
            Internal::rank(x, rx, TieMethod::TieAverage);
            Internal::rank(y, ry, TieMethod::TieAverage);
  
            return corrPearson(rx, ry);
        }
    }
}

#endif
