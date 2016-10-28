#include <catch.hpp>
#include "parsers/parser_exp.hpp"

using namespace Anaquin;

TEST_CASE("ParserExp_1")
{
    auto r = ParserExp::parse(Reader("tests/data/experiment.txt"));
    
    REQUIRE(r.samps.size() == 2);
    REQUIRE(r.samps[Mix_1].size() == 3);

    REQUIRE(r.samps[Mix_1][0].p1 == "L.R.1.1_val_1.fq");
    REQUIRE(r.samps[Mix_1][0].p2 == "L.R.1.2_val_2.fq");
    REQUIRE(r.samps[Mix_1][1].p1 == "L.R.2.1_val_1.fq");
    REQUIRE(r.samps[Mix_1][1].p2 == "L.R.2.2_val_2.fq");
    REQUIRE(r.samps[Mix_1][2].p1 == "L.R.3.1_val_1.fq");
    REQUIRE(r.samps[Mix_1][2].p2 == "L.R.3.2_val_2.fq");

    REQUIRE(r.samps[Mix_2][0].p1 == "L.R.4.1_val_1.fq");
    REQUIRE(r.samps[Mix_2][0].p2 == "L.R.4.2_val_2.fq");
    REQUIRE(r.samps[Mix_2][1].p1 == "L.R.5.1_val_1.fq");
    REQUIRE(r.samps[Mix_2][1].p2 == "L.R.5.2_val_2.fq");
    REQUIRE(r.samps[Mix_2][2].p1 == "L.R.6.1_val_1.fq");
    REQUIRE(r.samps[Mix_2][2].p2 == "L.R.6.2_val_2.fq");
}
