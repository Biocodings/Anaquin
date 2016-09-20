#include <catch.hpp>
#include "parsers/parser_exp.hpp"

using namespace Anaquin;

TEST_CASE("ParserExp_1")
{
    auto r = ParserExp::parse(Reader("tests/data/experiment.txt"));
    
    REQUIRE(r.samps.size() == 2);
    REQUIRE(r.samps[Mix_1].size() == 3);

    REQUIRE(r.samps[Mix_1][0].first  == "L.R.1.1_val_1.fq");
    REQUIRE(r.samps[Mix_1][0].second == "L.R.1.2_val_2.fq");
    REQUIRE(r.samps[Mix_1][1].first  == "L.R.2.1_val_1.fq");
    REQUIRE(r.samps[Mix_1][1].second == "L.R.2.2_val_2.fq");
    REQUIRE(r.samps[Mix_1][2].first  == "L.R.3.1_val_1.fq");
    REQUIRE(r.samps[Mix_1][2].second == "L.R.3.2_val_2.fq");
    
    REQUIRE(r.samps[Mix_2][0].first  == "L.R.4.1_val_1.fq");
    REQUIRE(r.samps[Mix_2][0].second == "L.R.4.2_val_2.fq");
    REQUIRE(r.samps[Mix_2][1].first  == "L.R.5.1_val_1.fq");
    REQUIRE(r.samps[Mix_2][1].second == "L.R.5.2_val_2.fq");
    REQUIRE(r.samps[Mix_2][2].first  == "L.R.6.1_val_1.fq");
    REQUIRE(r.samps[Mix_2][2].second == "L.R.6.2_val_2.fq");
}