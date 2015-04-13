#ifndef GI_PARSER_VCF_HPP
#define GI_PARSER_VCF_HPP

#include <set>
#include <string>
#include <vector>
#include <functional>
#include "types.hpp"

namespace Spike
{
    enum VCFVariantInfo
    {
        AC, // Allele frequency
    };
    
    enum Zygosity
    {
        Homozygous,
        Heterzygous,
    };

    typedef std::string Sequence;

    struct VCFVariant
    {
        // An identifier from the reference genome
        ChromoID chID;
        
        // The reference position, with the 1st base having position 1.
        BasePair pos;
        
        // Semi-colon separated list of unique identifiers where available
        VariantID varID;
        
        // Each base must be one of A,C,G,T,N (case insensitive). Multiple bases are permitted.
        Sequence ref;
        
        // Reference base - alternate non-reference alleles called on at least one of the samples.
        std::vector<Sequence> alts;
        
        Zygosity zy;
    };
    
    typedef std::function<void (const VCFVariant &)> VCFVariantF;
    
    struct ParserVCF
    {
        static void parse(const std::string &file, VCFVariantF);
    };
}

#endif