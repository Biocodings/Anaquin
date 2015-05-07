#include <catch.hpp>
#include "parsers/parser_vcf.hpp"

using namespace Spike;

TEST_CASE("Standard_DNA_1")
{
    std::vector<VCFVariant> vs;
    
    ParserVCF::parse("data/dna/variant.ChrT51.vcf", [&](const VCFVariant &v, const ParserProgress &)
    {
        vs.push_back(v);
    });

    REQUIRE(vs.size() == 245);
    
    REQUIRE(vs[0].l == Locus(373892, 373892));
    REQUIRE(vs[0].ref == "T");
    REQUIRE(vs[0].alt == "C");
    REQUIRE(vs[0].gt == HomozygousAlt);
    
    REQUIRE(vs[1].l == Locus(373899, 373899));
    REQUIRE(vs[1].ref == "TACTAACACGACGTGC");
    REQUIRE(vs[1].alt == "T");
    REQUIRE(vs[1].gt == HomozygousAlt);
}

TEST_CASE("Standard_DNA_2")
{
    std::vector<VCFVariant> vs;
    
    ParserVCF::parse("data/dna/hetero.ChrT51.vcf", [&](const VCFVariant &v, const ParserProgress &)
                     {
                         vs.push_back(v);
                     });

    REQUIRE(vs.size() == 139);
}
