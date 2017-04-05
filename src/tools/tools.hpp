#ifndef TOOLS_HPP
#define TOOLS_HPP

#include <map>
#include <string>
#include <memory>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <algorithm>

namespace Anaquin
{
    template <typename T> std::string toString(const T &x, unsigned n = 2)
    {
        std::ostringstream out;
        out << std::fixed << std::setprecision(n) << x;
        return out.str();
    }

    template <typename T> unsigned count(const std::map<T, unsigned> &x)
    {
        return std::accumulate(std::begin(x), std::end(x), 0, [](unsigned c, const std::pair<T, unsigned>& p)
        {
            return c + (p.second ? 1 : 0);
        });
    }
    
    template <typename X, typename F> unsigned countMap(const X &x, F f)
    {
        unsigned n = 0;
        
        for (const auto &i : x)
        {
            n += f(i.first, i.second);
        }
        
        return n;
    }
    
    template <typename T1, typename T2> T2 sum(const std::vector<T1> &x)
    {
        return std::accumulate(std::begin(x), std::end(x), T2(), [](T2 c, const T1 &p)
                               {
                                   return c + p;
                               });
    }
    
    template <typename T1, typename T2> T2 sum(const std::map<T1, T2> &x)
    {
        return std::accumulate(std::begin(x), std::end(x), T2(), [](T2 c, const std::pair<T1, T2>& p)
                               {
                                   return c + p.second;
                               });
    }
}

#endif
