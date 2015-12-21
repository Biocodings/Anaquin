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
    REQUIRE(r.n_chrT == 1200);

    REQUIRE(r.data.at(ChrT).eSN == 1.0);
    REQUIRE(r.data.at(ChrT).eSP == 1.0);
    REQUIRE(r.data.at(ChrT).bSN == 1.0);
    REQUIRE(r.data.at(ChrT).bSP == 1.0);
    REQUIRE(r.data.at(ChrT).tSN == 1.0);
    REQUIRE(r.data.at(ChrT).tSP == 1.0);
    REQUIRE(r.data.at(ChrT).iSN == Approx(0.996031746));
    REQUIRE(r.data.at(ChrT).iSP == Approx(0.996031746));
}