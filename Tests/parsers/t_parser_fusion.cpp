#include <catch.hpp>
#include "parsers/parser_fusion.hpp"

using namespace Anaquin;

TEST_CASE("Test_ParserFusion")
{
    std::vector<ParserFusion::Fusion> fs;
    
    ParserFusion::parse(Reader("tests/data/fusion/10K/fusions.out"), [&](const ParserFusion::Fusion &f, const ParserProgress &)
    {
        fs.push_back(f);
    });

    REQUIRE(fs[0].reads == 1);
    REQUIRE(fs[0].chr_1 == "chr1");
    REQUIRE(fs[0].chr_2 == "chr1");
    REQUIRE(fs[0].start_1 == 1217339);
    REQUIRE(fs[0].start_2 == 1217424);

    REQUIRE(fs[3].reads == 2);
    REQUIRE(fs[3].chr_1 == "chr1");
    REQUIRE(fs[3].chr_2 == "chr1");
    REQUIRE(fs[3].start_1 == 1405809);
    REQUIRE(fs[3].start_2 == 1489204);
}