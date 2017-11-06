#ifndef SS_STATS_HPP
#define SS_STATS_HPP

#include <map>
#include <vector>
#include <numeric>
#include <limits.h>
#include <algorithm>
#include <functional>
#include <ss/data/data.hpp>
#include <ss/data/errors.hpp>
#include <ss/internal/stats.hpp>

namespace SS
{
    template <typename T> typename T::value_type quant(const T &x, P probs)
    {
        return Internal::quant(x, probs);
    }
    
    template <typename T> typename T::value_type cov(const T &x, const T &y)
    {
        return Internal::cov(x, y);
    }

    template <class T> typename T::value_type var(const T &x)
    {
        return Internal::cov(x, x);
    }
    
    template <typename T> typename T::value_type SD(const T &x)
    {
        return sqrt(var(x));
    }
    
    template <typename T> typename T::value_type mean(const T &x)
    {
        return (1.0 * Internal::sum(x)) / Internal::count(x);
    }

    template <typename T> typename T::value_type med(T &x)
    {
        const auto n = x.size() / 2;
        std::nth_element(x.begin(), x.begin()+n, x.end());
        return x.at(n);
    }

    template <typename T1, typename T2> T2 med_(const std::map<T1, T2> &x)
    {
        std::vector<T2> t;
        
        for (const auto &i : x)
        {
            t.push_back(i.second);
        }
        
        return SS::med(t);
    }

    template <typename T> typename T::value_type min(const T &x)
    {
        return *(std::min_element(x.begin(), x.end()));
    }
    
    template <typename T> typename T::value_type max(const T &x)
    {
        return *(std::max_element(x.begin(), x.end()));
    }

    template <typename T> Real pearson(const T &x, const T &y)
    {
        return cov(x, y) / (SD(x) * SD(y));
    }
    
    enum TestType
    {
        Less,
        Greater,
        TwoSided,
    };
    
    namespace Internal
    {
        template <typename T> Real critical(T t, P conf, TestType type)
        {
            const auto alpha = 1.0 - conf;
            
            if (type == TwoSided)
            {
                return t(1.0 - (0.5 * alpha));
            }
            else if (type == Greater)
            {
                return t(1.0 - alpha);
            }
            else
            {
                return -t(1.0 - alpha);
            }
        }
        
        inline P pval(Real x, P cdf, TestType type)
        {
            if (type == TwoSided)
            {
                return 2.0 * (cdf > 0.5 ? 1.0 - cdf : cdf);
            }
            else if (type == Greater)
            {
                return 1.0 - cdf;
            }
            else
            {
                return cdf;
            }
        }
    }
}

#endif
