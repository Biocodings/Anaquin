#include "gtest/gtest.h"
#include "parser_vcf.hpp"

TEST(ref_variant, ParserVCFTest)
{
    unsigned i = 0;
    unsigned j = 0;
    
    ParserVCF::parse("/Users/tedwong/Sources/QA/data/DNA/variant.ChrT51.vcf", [&](const VCFHeader &h)
    {
        i++;
    }, [&](const VCFVariant &v)
    {
        j++;
    });

    ASSERT_EQ(0, i);
    ASSERT_EQ(245, j);
}