#include <catch.hpp>
#include "unit/test.hpp"
#include "trans/t_assembly.hpp"

using namespace Anaquin;

TEST_CASE("TAssembly_CompareWithItself")
{
    Test::transA();
    
    TAssembly::Options o;

    o.rChrT = "data/TransQuin/ATR001.v032.gtf";

    const auto r = TAssembly::analyze(o.rChrT, o);

    REQUIRE(r.chrT_exons == 1200);
    REQUIRE(r.endo_exons == 1200);
    REQUIRE(r.chrT_trans == 0);
    REQUIRE(r.endo_trans == 0);

    REQUIRE(r.data.at(ChrT).eSN == 1.0);
    REQUIRE(r.data.at(ChrT).eSP == 1.0);
    REQUIRE(r.data.at(ChrT).bSN == 1.0);
    REQUIRE(r.data.at(ChrT).bSP == 1.0);
    REQUIRE(r.data.at(ChrT).tSN == 1.0);
    REQUIRE(r.data.at(ChrT).tSP == 1.0);
    REQUIRE(r.data.at(ChrT).iSN == Approx(0.996031746));
    REQUIRE(r.data.at(ChrT).iSP == Approx(0.996031746));
}