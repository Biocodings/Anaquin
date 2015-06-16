//#include <catch.hpp>
//#include "parsers/parser_vcf.hpp"
//
//using namespace Spike;
//
//TEST_CASE("ParserVCF_Sample_File")
//{
//    std::vector<VCFVariant> vs;
//
//    ParserVCF::parse(Reader("data/dna/DNA.var.vcf"), [&](const VCFVariant &v, const ParserProgress &)
//    {
//        vs.push_back(v);
//    });
//
//    REQUIRE(vs.size() == 245);
//    
//    REQUIRE(vs[0].l  == Locus(373892, 373892));
//    REQUIRE(vs[0].r  == "T");
//    REQUIRE(vs[0].a  == "C");
//    REQUIRE(vs[0].gt == HomozygousAlt);
//    
//    REQUIRE(vs[1].l  == Locus(373899, 373899));
//    REQUIRE(vs[1].r  == "TACTAACACGACGTGC");
//    REQUIRE(vs[1].a  == "T");
//    REQUIRE(vs[1].gt == HomozygousAlt);
//}