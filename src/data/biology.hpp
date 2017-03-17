#ifndef BIOLOGY_HPP
#define BIOLOGY_HPP

#include <map>
#include <string>
#include "tools/errors.hpp"

namespace Anaquin
{
    /*
     * Complement DNA string
     */
    
    template <typename T> void complement(T &x)
    {
        std::map<char, char> m = { { 'A', 'T' },
                                   { 'T', 'A' },
                                   { 'G', 'C' },
                                   { 'C', 'G' },
                                   { 'N', 'N' } };

        std::transform(x.begin(), x.end(), x.begin(), [&](char c)
        {
            A_CHECK(m.count(c), "Unknown DNA base: " + std::to_string(c));
            return m[c];
        });
    }

    /*
     * Standardize chromosome name. For example, "1" and "chr1" should have the same meaning.
     */
    
    template <typename T> T standChr(const T &x)
    {
        return !x.empty() && isdigit(x.front()) ? "chr" + x : x;
    }
}

#endif
