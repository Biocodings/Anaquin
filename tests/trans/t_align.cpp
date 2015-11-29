#include <catch.hpp>
#include "unit/test.hpp"
#include "trans/t_align.hpp"

using namespace Anaquin;

TEST_CASE("TAlign_T_1000_First_100")
{
    Test::clear();
    
    const auto r1 = Test::test("-t TransAlign -m data/trans/MTR002.v013.csv -rgtf data/trans/ATR001.v032.gtf -usam tests/data/T_1000/B/accepted_hits_1000.sam");

    REQUIRE(r1.status == 0);

    Test::transA();
    const auto r2 = TAlign::stats("tests/data/T_1000/B/accepted_hits_1000.sam");

    REQUIRE(r2.unmapped  == 0);
    REQUIRE(r2.pe.m.sn() == Approx(0.9301025163));
    REQUIRE(r2.pe.m.sp() == Approx(0.9900793651));
    REQUIRE(r2.pi.m.sn() == Approx(0.8613678373));
    REQUIRE(r2.pi.m.sp() == Approx(1.0));
    REQUIRE(r2.pb.m.sn() == Approx(0.0012799979));
    REQUIRE(r2.pb.m.sp() == Approx(0.9845360825));
}