#ifndef GI_PARSER_VCF_HPP
#define GI_PARSER_VCF_HPP

#include <set>
#include <functional>
#include "types.hpp"

enum VCFVariantInfo
{
};

struct VCFVariant
{
    // An identifier from the reference genome
    ChromoID chromoID;
    
    // The reference position, with the 1st base having position 1.
    BasePair pos;
    
    // Semi-colon separated list of unique identifiers where available
    VariantID varID;
    
    // Each base must be one of A,C,G,T,N (case insensitive). Multiple bases are permitted.
    Sequence ref;
    
    // Reference base - alternate non-reference alleles called on at least one of the samples.
    std::set<Sequence> alts;
    
    // Whether this variant has been filtered
    bool filtered;
};

struct VCFHeader
{
    
};

typedef std::function<void (const VCFHeader &)> VCFHeaderF;
typedef std::function<void (const VCFVariant &)> VCFVariantF;

struct ParserVCF
{
    static void parse(const std::string &file, VCFHeaderF &, VCFVariantF &);
};

#endif