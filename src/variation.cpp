#include "variation.hpp"
#include "parser_vcf.hpp"
#include "standard_factory.hpp"

VariationStats Variation::analyze(const std::string &file)
{
    const auto r = StandardFactory::reference();

    ParserVCF::parse("/Users/tedwong/Sources/QA/data/DNA/variant.ChrT51.vcf", [&](const VCFVariant &v)
                     {
                         
                     });

    
    
    VariationStats stats;
    
    return stats;
}