#ifndef GI_DENOVO_HPP
#define GI_DENOVO_HPP

#include <vector>
#include "types.hpp"

namespace Spike
{
    struct Contig
    {
        BasePair l;
        std::string id;
        std::string seq;
    };
    
    struct DVRefStats
    {
        /*
         * The total number of aligned bases in the reference, divided by the genome size. A base in the reference genome
         * is counted as aligned if at least one contig has at least one alignment to this base. Contigs from repeated
         * regions may map to multiple places, and thus may be counted multiple times in this quantity.
         */

        double covered;

        /*
         * The total number of aligned bases in the assembly, divided by the total number of aligned bases in the reference.
         * If the assembly contains many contigs that cover the same regions of the reference, its duplication ratio may
         * be much > 1.
         */

        double duplicate;

        Counts aligned, unaligned;
    };

    struct DNStats
    {
        std::vector<Contig> contigs;

        BasePair min, max;
        BasePair mean, sum;
        BasePair N20, N50, N80;
    };

    struct DNAsssembly
    {
        static DNStats stats(const std::string &file);
    };
}

#endif