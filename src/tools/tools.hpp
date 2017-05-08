#ifndef TOOLS_HPP
#define TOOLS_HPP

#include <map>
#include <cmath>
#include <string>
#include <memory>
#include <vector>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/replace.hpp>

namespace Anaquin
{
    typedef std::string Tok;
    typedef std::vector<Tok> Toks;
    
    template <typename T> static void split(const Tok &x, const Tok &d, T &r)
    {
        r.clear();
        boost::split(r, x, boost::is_any_of(d));
    }

    inline Tok join(const Toks &x, const std::string &d)
    {
        return boost::algorithm::join(x, d);
    }
    
    // Eg: "C_12_D" to "C_12"
    inline Tok first(const Tok &x, const Tok &d)
    {
        Toks toks;
        toks.clear();
        split(x, d, toks);
        return toks.front();
    }

    // Eg: "C_12_D" to "D"
    inline Tok last(const Tok &x, const Tok &d)
    {
        Toks toks;
        toks.clear();
        split(x, d, toks);
        return toks.back();
    }
    
    // Eg: "C_12_D" to "C_12"
    inline Tok noLast(const Tok &x, const Tok &d)
    {
        Toks toks;
        split(x, d, toks);
        toks.pop_back();
        return join(toks, d);
    }

    inline std::string remove(const std::string &s1, const std::string &s2)
    {
        auto x = s1;
        boost::replace_all(x, s2, "");
        return x;
    }

    inline bool hasSub(const std::string &s1, const std::string &s2)
    {
        return s1.find(s2) != std::string::npos;
    }

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
    
    inline long double ss2ld(const std::string &x)
    {
        if (x == "NA" || x == "-" || x == "*")
        {
            return NAN;
        }
        
        std::istringstream os(x);
        long double p;
        os >> p;
        return p;
    }
    
    // Convert floating number to scientific notation
    inline std::string ld2ss(long double p)
    {
        if (isnan(p) || !isfinite(p))
        {
            return "-";
        }

        std::ostringstream out;
        
        // Convert the probability into scientific notation
        out << std::scientific << p;
        
        return out.str();
    }
    
    template <typename T> std::string toString(const T &x, unsigned n = 2)
    {
        if (isnan(x) || !isfinite(x))
        {
            return "-";
        }
        
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
