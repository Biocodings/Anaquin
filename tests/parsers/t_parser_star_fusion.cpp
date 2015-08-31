#include <catch.hpp>
#include "parsers/parser_star_fusion.hpp"

using namespace Anaquin;

TEST_CASE("ParserStarFusion_Test")
{
    std::vector<ParserStarFusion::Fusion> fs;

    ParserStarFusion::parse(Reader("tests/data/F_1000/star-fusion.fusion_candidates.txt"), [&](const ParserStarFusion::Fusion &f, const ParserProgress &)
    {
        fs.push_back(f);
    });

    REQUIRE(fs.size() == 638);
}