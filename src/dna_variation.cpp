#include "parser_vcf.hpp"
#include "dna_variation.hpp"
#include "standard_factory.hpp"

VariationStats DNAVariation::analyze(const std::string &file)
{
    const auto r = StandardFactory::reference();

    ParserVCF::parse("data/DNA/variant.ChrT51.vcf", [&](const VCFVariant &v)
                     {
                         
                     });

    
    
    VariationStats stats;
    
    return stats;
}