#ifndef CONVERT_HPP
#define CONVERT_HPP

#include <string>

namespace Anaquin
{
    inline long double ns2ld(const std::string &s)
    {
        std::istringstream os(s);
        long double p;
        os >> p;
        return p;
    }
    
    inline std::string ld2ns(long double p)
    {
        std::ostringstream out;
        
        // Convert the probability into scientific notation
        out << std::scientific << p;
        
        return out.str();
    }
}

#endif