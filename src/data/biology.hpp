#ifndef BIOLOGY_HPP
#define BIOLOGY_HPP

#include <map>
#include <string>

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

    /*
     * Standardize chromosome names. For example, "1" and "chr1" should mean the same thing.
     */
    
    template <typename T> T standChr(const T &cID)
    {
        return cID.find("chr") == std::string::npos ? "chr" + cID : cID;
    }
}

#endif
