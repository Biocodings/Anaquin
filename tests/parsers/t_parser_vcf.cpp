#include <catch.hpp>
#include "data/reader.hpp"
#include "parsers/parser_vcf.hpp"

using namespace Anaquin;

//TEST_CASE("ParserVCF_Invalid")
//{
//    REQUIRE(!ParserVCF::isVCF(Reader("tests/data/Invalid.vcf")));
//}

TEST_CASE("ParserVCF_AVA026")
{
    std::vector<ParserVCF::Data> x;
    
    ParserVCF::parse(Reader("tests/data/AVA026_v001.vcf"), [&](const ParserVCF::Data &d, const ParserProgress &)
    {
        x.push_back(d);
    });

    REQUIRE(x.size() == 102);

    REQUIRE(x[0].cID == "chr1");
    REQUIRE(x[0].l.start == 58701656);
    REQUIRE(x[0].id  == "GI_086");
    REQUIRE(x[0].ref == "GGAAT");
    REQUIRE(x[0].alt == "G");
    REQUIRE(x[0].type() == Mutation::Deletion);
    
    REQUIRE(x[1].cID == "chr1");
    REQUIRE(x[1].l.start == 60005804);
    REQUIRE(x[1].id  == "GI_048");
    REQUIRE(x[1].ref == "C");
    REQUIRE(x[1].alt == "CTT");
    REQUIRE(x[1].type() == Mutation::Insertion);
}
