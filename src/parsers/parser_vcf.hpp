#ifndef GI_PARSER_VCF_HPP
#define GI_PARSER_VCF_HPP

#include <set>
#include "standard.hpp"
#include "parsers/parser.hpp"

namespace Spike
{
    enum VCFVariantInfo
    {
        AC, // Allele frequency
    };

    struct VCFVariant
    {
        ChromoID id;

        // The reference position, with the 1st base having position 1
        Locus l;

        Mutation m;

        Sequence ref;
        Sequence alt;

        Genotype gt;
        
        // Allelle frequency
        Counts af;

        // Allele count in genotypes
        Counts ac;

        // Total number of alleles in called genotypes
        Counts an;

        // Combined depth across samples
        unsigned dp;
    };

    struct ParserVCF
    {
        static void parse(const std::string &file, std::function<void (const VCFVariant &, const ParserProgress &)>);
    };
}

#endif