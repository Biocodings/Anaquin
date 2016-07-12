#ifndef CONVERT_HPP
#define CONVERT_HPP

#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>

namespace Anaquin
{
    inline long double ss2ld(const std::string &s)
    {
        std::istringstream os(s);
        long double p;
        os >> p;
        return p;
    }
    
    inline std::string ld2ss(long double p)
    {
        std::ostringstream out;
        
        // Convert the probability into scientific notation
        out << std::scientific << p;
        
        return out.str();
    }

    template <typename T> std::string x2ns(const T &x)
    {
        auto f = [&](unsigned n=2)
        {
            std::ostringstream out;
            out << std::fixed << std::setprecision(n) << x;
            return out.str();
        };

        return isnan(x) ? "NA" : f(x);
    }
}

#endif