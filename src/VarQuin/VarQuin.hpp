#ifndef VARQUIN_HPP
#define VARQUIN_HPP

#include "data/data.hpp"
#include "data/standard.hpp"
#include "parsers/parser_vcf.hpp"
#include "parsers/parser_variants.hpp"
#include <boost/algorithm/string/predicate.hpp>

namespace Anaquin
{
    struct VariantStats
    {
        // Number of SNPs detected
        Counts n_snp;

        // Number of indels detected
        Counts n_ind;
    };
    
    struct VariantMatch
    {
        // The called variant
        CalledVariant query;

        // Matched by position?
        const Variant *match = nullptr;
        
        // Matched by variant allele? Only if position is matched.
        bool alt;
        
        // Matched by reference allele? Only if position is matched.
        bool ref;
    };

    // Eg: chrev1, chrev10 etc...
    inline bool isReverseGenome(const ChrID &cID)
    {
        A_ASSERT(!cID.empty());
        return cID.find("rev") != std::string::npos;
    }

    inline ChrID toReverse(const ChrID &cID)
    {
        auto x = cID;
        
        // Eg; chr2 to chrev2
        boost::replace_all(x, "chr", "chrev");
        
        A_ASSERT(!x.empty());
        
        return x;
    }
}

#endif
