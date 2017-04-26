#include <catch.hpp>
#include "parsers/parser_salmon.hpp"

using namespace Anaquin;

TEST_CASE("ParserSalmon_1")
{
    std::vector<ParserSalmon::Data> x;
    
    ParserSalmon::parse(Reader("tests/data/quant.sf"), [&](const ParserSalmon::Data &d, const ParserProgress &)
    {
        x.push_back(d);
    });

    REQUIRE(x.size() == 200);

    REQUIRE(x[0].name  == "CI_012_R");
    REQUIRE(x[0].abund == 18.5236);

    REQUIRE(x[1].name  == "CI_015_R");
    REQUIRE(x[1].abund == 12.5155);

    REQUIRE(x[2].name  == "CI_019_R");
    REQUIRE(x[2].abund == 37.582);
}
