//#include <catch.hpp>
//#include "unit/test.hpp"
//#include "variant/v_align.hpp"
//
//using namespace Anaquin;
//
//TEST_CASE("VAlign_GM_VARMXA")
//{
//    Test::variant();
//
//    const auto r = VAlign::analyze("tests/data/GM_VARMXA_CONA/aligned.sam");
//
//    REQUIRE(r.p.m.sp() == Approx(0.9463336876));
//    REQUIRE(r.p.m.sn() == Approx(0.961663067));
//    REQUIRE(r.p.m.nq == 1882);
//    REQUIRE(r.p.m.nr == 1852);
//    REQUIRE(r.p.s.id == "D_1_3_V");
//    REQUIRE(r.p.s.counts == 1781);
//    REQUIRE(r.p.s.abund == Approx(20000));
//}