#ifndef BIOLOGY_HPP
#define BIOLOGY_HPP

#include <map>

namespace Anaquin
{
    template <typename T> void complement(T &str)
    {
        std::map<char, char> m = { { 'A', 'T' },
                                   { 'T', 'A' },
                                   { 'G', 'C' },
                                   { 'C', 'G' },
                                   { 'N', 'N' } };
        for (auto &i : str)
        {
            if (!m.count(i))
            {
                throw std::runtime_error("Unknown DNA base: " + str);
            }
            
            i = m[i];
        }
    }
}

#endif
