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

    REQUIRE(r.chrT->n_expT == 0);
    REQUIRE(r.chrT->n_chrT == 2730);

    /*
     * It's very likely cuffcompare is just wrong. Everything here should be 1.0.
     */
    
    REQUIRE(r.chrT->exonSN   == 1.0);
    REQUIRE(r.chrT->exonSP   == 1.0);
    REQUIRE(r.chrT->intronSN == Approx(0.996031746));
    REQUIRE(r.chrT->intronSP == Approx(0.996031746));
    REQUIRE(r.chrT->baseSN   == 1.0);
    REQUIRE(r.chrT->baseSP   == 1.0);
    REQUIRE(r.chrT->transSN  == 1.0);
    REQUIRE(r.chrT->transSP  == 1.0);
}