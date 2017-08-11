#ifndef TOOLS_HPP
#define TOOLS_HPP

#include <set>
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
#include <boost/algorithm/string/predicate.hpp>

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
    
    inline bool isEnded(const Tok &x, const Tok &y)
    {
        return boost::algorithm::ends_with(x, y);
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
        if (std::isnan(p) || !std::isfinite(p))
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
        if (std::isnan(x) || !std::isfinite(x))
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
    
    template <typename T> double mean(const T &x)
    {
        return std::accumulate(x.begin(), x.end(), 0.0) / x.size();
    }

    template <typename T1, typename T2> T2 sum(const std::vector<T1> &x)
    {
        return std::accumulate(std::begin(x), std::end(x), T2(), [](T2 c, const T1 &p)
        {
            return c + p;
        });
    }

    template <typename T1, typename T2> unsigned nonZero(const std::map<T1, T2> &x)
    {
        unsigned n = 0;

        for (const auto &i : x)
        {
            if (i.second)
            {
                n++;
            }
        }
        
        return n;
    }
    
    typedef std::string Path;
    typedef std::string FileName;
    
    inline FileName path2file(const Path &path)
    {
        auto tmp = path;
        const auto last = path.find_last_of("\\/");
        
        if (std::string::npos != last)
        {
            tmp.erase(0, last + 1);
        }
        
        return tmp;
    }
    
    inline std::vector<FileName> path2file(const std::vector<FileName> &files)
    {
        auto tmp = files;
        
        for (auto i = 0; i < tmp.size(); i++)
        {
            tmp[i] = path2file(tmp[i]);
        }
        
        return tmp;
    }

    template <typename T1, typename T2> T2 sum(const std::map<T1, T2> &x)
    {
        return std::accumulate(std::begin(x), std::end(x), T2(), [](T2 c, const std::pair<T1, T2>& p)
        {
            return c + p.second;
        });
    }

    template <typename Key, typename Value> std::set<Key> keys(const std::map<Key, Value> &x)
    {
        std::set<Key> keys;
        
        for (const auto i: x)
        {
            keys.insert(i.first);
        }
        
        return keys;
    }
}

#endif
