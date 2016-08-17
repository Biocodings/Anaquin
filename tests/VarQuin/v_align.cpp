//#include <catch.hpp>
//#include "test.hpp"
//#include "VarQuin/v_align.hpp"
//
//using namespace Anaquin;
//
//TEST_CASE("VAlign_Test1")
//{
//    Test::VarQuinBed();
//    const auto r = VAlign::analyze("tests/data/test1.bam");
//
//    REQUIRE(r.countSyn  == 356332);
//    REQUIRE(r.countGen  == 0);
//    REQUIRE(r.n_unmap == 2015);
//    
////    REQUIRE(r.p.m.ac() == Approx(0.9463336876));
////    REQUIRE(r.p.m.sn() == Approx(0.961663067));
////    REQUIRE(r.p.m.nq == 1882);
////    REQUIRE(r.p.m.nr == 1852);
////    REQUIRE(r.p.s.id == "D_1_3_V");
////    REQUIRE(r.p.s.counts == 1781);
////    REQUIRE(r.p.s.abund == Approx(20000));
//}