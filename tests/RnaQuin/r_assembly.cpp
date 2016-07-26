#include <catch.hpp>
#include "unit/test.hpp"
#include "RnaQuin/r_assembly.hpp"

using namespace Anaquin;

// Defined in main.cpp
extern void SetGTFRef(const FileName &);

TEST_CASE("RAssembly_CompareWithItself")
{
    Test::transA();
    
    const auto file = "data/RnaQuin/ARN020_v001.gtf";

    SetGTFRef(file);
    const auto r = RAssembly::analyze(file);
    SetGTFRef("");

    REQUIRE(r.gExons == 0);
    REQUIRE(r.gGenes == 0);
    REQUIRE(r.gTrans == 0);
    REQUIRE(r.gIntrs == 0);
    REQUIRE(r.sExons == 869);
    REQUIRE(r.sGenes == 754);
    REQUIRE(r.sTrans == 164);
    REQUIRE(r.sIntrs == 78);

    REQUIRE(r.data.at(ChrT).eSN == 1.0);
    REQUIRE(r.data.at(ChrT).eSP == 1.0);
    REQUIRE(r.data.at(ChrT).bSN == 1.0);
    REQUIRE(r.data.at(ChrT).bSP == 1.0);
    REQUIRE(r.data.at(ChrT).tSN == 1.0);
    REQUIRE(r.data.at(ChrT).tSP == 1.0);
    REQUIRE(r.data.at(ChrT).iSN == Approx(1.0));
    REQUIRE(r.data.at(ChrT).iSP == Approx(1.0));
}