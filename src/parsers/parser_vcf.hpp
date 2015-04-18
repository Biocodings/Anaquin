#ifndef GI_PARSER_VCF_HPP
#define GI_PARSER_VCF_HPP

#include <set>
#include "types.hpp"
#include "standard.hpp"

namespace Spike
{
    enum VCFVariantInfo
    {
        AC, // Allele frequency
    };

    struct VCFVariant
    {
        // An identifier from the reference genome
        ChromoID chID;

        // The reference position, with the 1st base having position 1
        BasePair pos;

        // Semi-colon separated list of unique identifiers where available
        VariantID varID;
        
        Sequence r;
        
        // Reference base - alternate non-reference alleles called on at least one of the samples
        std::set<Sequence> alts;

        Zygosity zy;
    };

    typedef std::function<void (const VCFVariant &)> VCFVariantF;

    struct ParserVCF
    {
        static void parse(const std::string &file, VCFVariantF);
    };
}

#endif