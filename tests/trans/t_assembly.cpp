#include <catch.hpp>
#include "unit/test.hpp"
#include "trans/t_assembly.hpp"

using namespace Anaquin;

TEST_CASE("TAssembly_Reference")
{
    Test::transA();
    
    TAssembly::Options o;
    
    o.ref   = "data/trans/ATR001.v032.gtf";
    o.query = "data/trans/ATR001.v032.gtf";

    const auto r = TAssembly::analyze("data/trans/ATR001.v032.gtf", o);

    REQUIRE(r.n_expT == 0);
    REQUIRE(r.n_chrT == 2730);

    /*
     * It's very likely cuffcompare is just wrong. Everything here should be 1.0.
     */
    
    REQUIRE(r.exonSN   == 1.0);
    REQUIRE(r.exonSP   == 1.0);
    REQUIRE(r.intronSN == Approx(0.996031746));
    REQUIRE(r.intronSP == Approx(0.996031746));
    REQUIRE(r.baseSN   == 1.0);
    REQUIRE(r.baseSP   == 1.0);
    REQUIRE(r.transSN  == 1.0);
    REQUIRE(r.transSP  == 1.0);
}