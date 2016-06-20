#ifndef BIOLOGY_HPP
#define BIOLOGY_HPP

#include <map>
#include "data/types.hpp"

namespace Anaquin
{
    enum Strand
    {
        Forward,
        Backward,
    };

    enum RNAFeature
    {
        Exon,
        Gene,
        Intron,
        Transcript,
    };
    
    enum Mutation
    {
        SNP,
        Insertion,
        Deletion
    };
    
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