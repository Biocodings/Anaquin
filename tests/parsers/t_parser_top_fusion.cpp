#include <catch.hpp>
#include "parsers/parser_top_fusion.hpp"

using namespace Anaquin;

TEST_CASE("Test_ParserFusion")
{
    std::vector<ParserTopFusion::Fusion> fs;

    ParserTopFusion::parse(Reader("tests/examples/fusions.out"), [&](const ParserTopFusion::Fusion &f, const ParserProgress &)
    {
        fs.push_back(f);
    });

    REQUIRE(fs[0].reads == 1);
    REQUIRE(fs[0].chr_1 == "chr1");
    REQUIRE(fs[0].chr_2 == "chr1");
    REQUIRE(fs[0].l1 == 1217339);
    REQUIRE(fs[0].l2 == 1217424);

    REQUIRE(fs[3].reads == 2);
    REQUIRE(fs[3].chr_1 == "chr1");
    REQUIRE(fs[3].chr_2 == "chr1");
    REQUIRE(fs[3].l1 == 1405809);
    REQUIRE(fs[3].l2 == 1489204);
}