#ifndef CONVERT_HPP
#define CONVERT_HPP

#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>

namespace Anaquin
{
    inline double s2d(const std::string &x)
    {
        return (x == "NA" || x == "-") ? NAN : stod(x);
    }

    inline long double ss2ld(const std::string &s)
    {
        if (s == "NA" || s == "-" || s == "*")
        {
            return NAN;
        }

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
        return isnan(x) ? "NA" : std::to_string(x);
    }
}

#endif