#include <catch.hpp>
#include "rna/r_assembly.hpp"

using namespace Anaquin;

TEST_CASE("RAssembly_Simulation")
{
    const auto r = RAssembly::analyze("tests/data/rna/A1/transcripts.gtf");

    REQUIRE(r.pe.m.nq == 1200);
    REQUIRE(r.pe.m.nr == 1199);
    REQUIRE(r.pe.m.sp() == Approx(0.9975));
    REQUIRE(r.pe.m.sn() == Approx(0.9983319433));
    REQUIRE(r.pe.s.id == "R2_38_1");
    REQUIRE(r.pe.s.counts == 1);
    REQUIRE(r.pe.s.abund == Approx(1.0));

    REQUIRE(r.pi.m.nq == 1035);
    REQUIRE(r.pi.m.nr == 887);
    REQUIRE(r.pi.m.sp() == Approx(0.8415458937));
    REQUIRE(r.pi.m.sn() == Approx(0.9819616685));
    REQUIRE(r.pi.s.id == "R2_38_4");
    REQUIRE(r.pi.s.counts == 1);
    REQUIRE(r.pi.s.abund == Approx(4.0));
    
    REQUIRE(r.pb.m.nq == 151192);
    REQUIRE(r.pb.m.nr == 149219);
    REQUIRE(r.pb.m.sp() == Approx(0.9869503677));
    REQUIRE(r.pb.m.sn() == Approx(1.0));
    REQUIRE(r.pb.s.id == "R2_28");
    REQUIRE(r.pb.s.counts == 1);
    REQUIRE(r.pb.s.abund == Approx(57.0));

    REQUIRE(r.pt.m.nq == 165);
    REQUIRE(r.pt.m.nr == 164);
    REQUIRE(r.pt.m.sp() == Approx(0.9818181818));
    REQUIRE(r.pt.m.sn() == Approx(0.987804878));
    REQUIRE(r.pt.s.id == "R2_76_1");
    REQUIRE(r.pt.s.counts == 1);
    REQUIRE(r.pt.s.abund == Approx(1.0));
}