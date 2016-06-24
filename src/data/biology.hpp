#ifndef BIOLOGY_HPP
#define BIOLOGY_HPP

#include <map>

namespace Anaquin
{
    template <typename T> void complement(T &s)
    {
        std::map<char, char> m = { { 'A', 'T' },
                                   { 'T', 'A' },
                                   { 'G', 'C' },
                                   { 'C', 'G' } };
        for (auto &i : s)
        {
            if (!m.count(i))
            {
                throw "Unknown DNA base: " + std::to_string(i);
            }
            
            i = m[i];
        }
    }
}

#endif