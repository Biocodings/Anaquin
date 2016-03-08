#include <catch.hpp>
#include "parsers/parser_star_fusion.hpp"

using namespace Anaquin;

TEST_CASE("ParserStarFusion_Test")
{
    std::vector<ParserStarFusion::Data> fs;

    ParserStarFusion::parse(Reader("tests/data/star-fusion.fusion_candidates.txt"), [&](const ParserStarFusion::Data &f, const ParserProgress &)
    {
        fs.push_back(f);
    });

    REQUIRE(fs.size() == 638);
}