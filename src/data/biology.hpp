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
}

#endif