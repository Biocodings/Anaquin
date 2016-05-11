#include <catch.hpp>
#include "data/reader.hpp"
#include "parsers/parser_stamp.hpp"

using namespace Anaquin;

TEST_CASE("ParserStamp_Test")
{
    std::vector<ParserStamp::Data> x;
    
    ParserStamp::parse(Reader("tests/data/STAMP.tsv"), [&](const ParserStamp::Data &d, const ParserProgress &)
    {
        x.push_back(d);
    });
    
    REQUIRE(x.size() == 3);
    
    REQUIRE(x[0].id == "S1");
    REQUIRE(x[1].id == "S2");
    REQUIRE(x[2].id == "S3");

    REQUIRE(x[0].p == Approx(0.135574438411));
    REQUIRE(x[1].p == Approx(0.349702005061));
    REQUIRE(x[2].p == Approx(0.219462971344));

    REQUIRE(x[0].q == Approx(0.135574438411));
    REQUIRE(x[1].q == Approx(0.349702005061));
    REQUIRE(x[2].q == Approx(0.219462971344));

    REQUIRE(x[0].effect == Approx(0.465179844804));
    REQUIRE(x[1].effect == Approx(0.218664369634));
    REQUIRE(x[2].effect == Approx(0.345966676945));
}