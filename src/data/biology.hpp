#ifndef BIOLOGY_HPP
#define BIOLOGY_HPP

#include <map>
#include "tools/errors.hpp"

namespace Anaquin
{
    enum Strand
    {
        Forward,
        Backward,
        Either,
    };
    
    enum class RNAFeature
    {
        Exon,
        Gene,
        Intron,
        Transcript,
    };

    enum Variation
    {
        SNP,
        Insertion,
        Deletion,
        Inversion,
        Duplication,
    };

    enum class Genotype
    {
        Somatic,
        Homozygous,
        Heterzygous
    };
    
    enum class Mutation
    {
        Germline,
        Somatic
    };
    
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
            A_ASSERT(m.count(c));
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
