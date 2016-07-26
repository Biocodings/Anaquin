#ifndef CONVERT_HPP
#define CONVERT_HPP

#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>
#include <stdexcept>

namespace Anaquin
{
    inline double s2d(const std::string &x)
    {
        if (x == "NA" || x == "-")
        {
            return NAN;
        }
        
        try
        {
            return stod(x);
        }
        catch(...)
        {
            throw std::runtime_error("Failed to parse \"" + x + "\". This is not a number.");
        }
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
        if (isnan(p) || !std::isfinite(p))
        {
            return "NA";
        }
        
        std::ostringstream out;
        
        // Convert the probability into scientific notation
        out << std::scientific << p;
        
        return out.str();
    }

    template <typename T> std::string x2ns(const T &x)
    {
        if (isnan(x) || !std::isfinite(x))
        {
            return "NA";
        }
        
        return std::to_string(x);
    }
}

#endif