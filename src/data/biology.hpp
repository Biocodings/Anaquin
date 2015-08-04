#ifndef GI_BIOLOGY_HPP
#define GI_BIOLOGY_HPP

namespace Anaquin
{
    enum Strand
    {
        Forward,
        Backward,
    };

    enum Orientation
    {
        ForwardForward,
        ForwardReverse,
        ReverseForward,
        ReverseReverse
    };

    struct Contig
    {
        std::string id;

        // The sequence being assembled
        std::string seq;

        // Coverage in k-mer
        double k_cov;
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

    struct FusionBreak
    {
        inline bool operator<(const FusionBreak &x)  const { return id < x.id;  }
        inline bool operator==(const FusionBreak &x) const { return id == x.id; }

        // Where this fusion belongs
        SequinID id;

        // The position of the break-point
        Locus l;

        Orientation orient;
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