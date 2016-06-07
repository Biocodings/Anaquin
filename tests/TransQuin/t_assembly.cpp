//#include <catch.hpp>
//#include "unit/test.hpp"
//#include "TransQuin/t_assembly.hpp"
//
//using namespace Anaquin;
//
//TEST_CASE("TAssembly_CompareWithItself")
//{
//    Test::transA();
//    
//    TAssembly::Options o;
//
//    o.rAnnot = "data/TransQuin/ATR001.v032.gtf";
//
//    const auto r = TAssembly::analyze(o.rAnnot, o);
//
//    REQUIRE(r.cExons == 1200);
//    REQUIRE(r.eExons == 1200);
//    REQUIRE(r.cTrans == 0);
//    REQUIRE(r.eTrans == 0);
//
//    REQUIRE(r.data.at(ChrT).eSN == 1.0);
//    REQUIRE(r.data.at(ChrT).eSP == 1.0);
//    REQUIRE(r.data.at(ChrT).bSN == 1.0);
//    REQUIRE(r.data.at(ChrT).bSP == 1.0);
//    REQUIRE(r.data.at(ChrT).tSN == 1.0);
//    REQUIRE(r.data.at(ChrT).tSP == 1.0);
//    REQUIRE(r.data.at(ChrT).iSN == Approx(0.996031746));
//    REQUIRE(r.data.at(ChrT).iSP == Approx(0.996031746));
//}