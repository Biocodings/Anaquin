#ifndef GI_BIOLOGY_HPP
#define GI_BIOLOGY_HPP

namespace Spike
{
    enum Strand
    {
        Forward,
        Backward,
    };

    struct Contig
    {
        std::string id;

        // The sequence being assembled
        std::string seq;

        // Coverage in k-mer
        double k_cov;
    };

    enum RNALevel
    {
        Gene,
        Isoform
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