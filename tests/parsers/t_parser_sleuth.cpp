#include <catch.hpp>
#include "test.hpp"
#include "data/reader.hpp"
#include "parsers/parser_sleuth.hpp"

using namespace Anaquin;

TEST_CASE("ParserSleuth_1")
{
    Test::transA();

    std::vector<ParserSleuth::Data> x;
    
    ParserSleuth::parse(Reader("tests/data/sleuth.csv"), [&](const ParserSleuth::Data &d, const ParserProgress &)
    {
        x.push_back(d);
    });

    REQUIRE(x.size() == 164);

    REQUIRE(x[0].cID    == "chrIS");
    REQUIRE(x[0].gID    == "R2_68");
    REQUIRE(x[0].iID    == "R2_68_2");
    REQUIRE(x[0].logF_  == Approx(-4.378098453));
    REQUIRE(x[0].logFSE == Approx(0.093369275));
    REQUIRE(x[0].mean   == Approx(1.495902046));
    REQUIRE(x[0].p      == 0.0);
    REQUIRE(x[0].q      == 0.0);
    REQUIRE(isnan(x[0].samp1));
    REQUIRE(isnan(x[0].samp2));
    REQUIRE(x[0].status == DiffTest::Status::Tested);

    REQUIRE(x[17].cID    == "chrIS");
    REQUIRE(x[17].gID    == "R2_115");
    REQUIRE(x[17].iID    == "R2_115_2");
    REQUIRE(x[17].logF_  == Approx(2.357638407));
    REQUIRE(x[17].logFSE == Approx(0.169674522));
    REQUIRE(x[17].mean   == Approx(6.987154919));
    REQUIRE(x[17].p      == Approx(6.78625e-44));
    REQUIRE(x[17].q      == Approx(3.95864400992919e-44));
    REQUIRE(isnan(x[17].samp1));
    REQUIRE(isnan(x[17].samp2));
    REQUIRE(x[17].status == DiffTest::Status::Tested);
}