#include <catch.hpp>
#include "unit/test.hpp"
#include "trans/t_assembly.hpp"

using namespace Anaquin;

TEST_CASE("TAssembly_Command")
{
    Test::trans();

    const auto r = Test::test("-t TransAssembly -m data/trans/RNA_4_1.csv -rgtf data/trans/RNA_1.gtf -ugtf tests/data/trans/A1/transcripts.gtf");

    REQUIRE(r.status == 0);
    REQUIRE(r.output.find("Transcriptome Analysis") != std::string::npos);
}

TEST_CASE("TAssembly_Test")
{
    Test::trans();

    const auto r = TAssembly::analyze("tests/data/trans/A1/transcripts.gtf");

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