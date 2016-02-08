#include <catch.hpp>
#include "parsers/parser_top_fusion.hpp"

using namespace Anaquin;

TEST_CASE("Test_ParserFusion")
{
    std::vector<ParserTopFusion::Data> data;

    ParserTopFusion::parse(Reader("tests/examples/fusions.out"), [&](const ParserTopFusion::Data &x, const ParserProgress &)
    {
        data.push_back(x);
    });

    REQUIRE(data[0].reads == 1);
    REQUIRE(data[0].cID_1 == "chr1");
    REQUIRE(data[0].cID_2 == "chr1");
    REQUIRE(data[0].l1 == 1217339);
    REQUIRE(data[0].l2 == 1217424);

    REQUIRE(data[3].reads == 2);
    REQUIRE(data[3].cID_1 == "chr1");
    REQUIRE(data[3].cID_2 == "chr1");
    REQUIRE(data[3].l1 == 1405809);
    REQUIRE(data[3].l2 == 1489204);
}