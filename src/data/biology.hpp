#ifndef BIOLOGY_HPP
#define BIOLOGY_HPP

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
        Gene,
        Intron,
        StopCodon,
        StartCodon,
        Transcript,
    };
    
    enum VarType
    {
        SNP,
        Insertion,
        Deletion
    };

    enum Genotype
    {
        Heterzygous,
        HomozygousR,
        HomozygousA,
    };
}

#endif