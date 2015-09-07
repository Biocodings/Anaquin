#ifndef GI_BIOLOGY_HPP
#define GI_BIOLOGY_HPP

#include "data/types.hpp"

namespace Anaquin
{
    enum Strand
    {
        Forward,
        Backward,
    };

    struct Contig
    {
        ContigID id;

        // The sequence being assembled
        Sequence seq;

        // Coverage in k-mer
        Coverage k_cov;
    };

    enum RNAFeature
    {
        CDS,
        Exon,
        Intron,
        StopCodon,
        StartCodon,
        Transcript,
    };
    
    enum Mutation
    {
        SNP,
        Insertion,
        Deletion
    };

    enum Genotype
    {
        Heterzygous,
        HomozygousRef,
        HomozygousAlt,
    };

    struct FusionPoint
    {
        inline bool operator<(const FusionPoint &x)  const { return id < x.id;  }
        inline bool operator==(const FusionPoint &x) const { return id == x.id; }

        // Where this fusion belongs
        SequinID id;

        // The position of the break-point
        Base l1, l2;

        // Orientation for each of the segment
        Strand s1, s2;
    };

    template <typename Iter, typename T, typename F> bool find(const Iter &begin, const Iter &end, const T &t, F &r)
    {
        for (auto i = begin; i < end; i++)
        {
            if (i->l.contains(t.l))
            {
                r = *i;
                return true;
            }
        }
        
        return false;
    }
}

#endif