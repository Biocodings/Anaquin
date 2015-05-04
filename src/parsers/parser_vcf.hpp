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
        // An identifier from the reference genome
        ChromoID id;

        // The reference position, with the 1st base having position 1
        Locus l;

        Sequence r;

        /*
         * List of alternate non-reference alleles called on at least one of the samples
         */

        std::set<Sequence> alts;

        /*
         * Additional information
         */
        
        std::map<std::string, std::string> info;
        
        Zygosity zy;
    };

    struct ParserVCF
    {
        static void parse(const std::string &file, std::function<void (const VCFVariant &, const ParserProgress &)>);
    };
}

#endif