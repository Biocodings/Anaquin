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
    REQUIRE(r.sGenes == 78);
    REQUIRE(r.sTrans == 164);
    REQUIRE(r.sIntrs == 754);

    REQUIRE(r.data.at(ChrIS).eSN == 1.0);
    REQUIRE(r.data.at(ChrIS).eSP == 1.0);
    REQUIRE(r.data.at(ChrIS).bSN == 1.0);
    REQUIRE(r.data.at(ChrIS).bSP == 1.0);
    REQUIRE(r.data.at(ChrIS).tSN == 1.0);
    REQUIRE(r.data.at(ChrIS).tSP == Approx(0.905982906));
    REQUIRE(r.data.at(ChrIS).iSN == Approx(0.9960212202));
    REQUIRE(r.data.at(ChrIS).iSP == Approx(0.9960212202));
}