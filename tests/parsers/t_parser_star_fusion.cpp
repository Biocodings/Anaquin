#include <catch.hpp>
#include "parsers/parser_star_fusion.hpp"

using namespace Anaquin;

TEST_CASE("ParserSFusion_100K")
{
    std::vector<ParserStarFusion::Fusion> fs;

    ParserStarFusion::parse(Reader("tests/data/fusion/100K/star-fusion.fusion_candidates.txt"), [&](const ParserStarFusion::Fusion &f, const ParserProgress &)
    {
        fs.push_back(f);
    });

    REQUIRE(fs.size() == 638);
}