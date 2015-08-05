#include <catch.hpp>
#include "unit/test.hpp"
#include "trans/t_align.hpp"

using namespace Anaquin;

TEST_CASE("TAlign_Cufflinks")
{
    Test::trans();

    // The sample file was taken from Cufflink's source distribution
    const auto r = TAlign::analyze("tests/data/trans/cufflinks.sam");

    REQUIRE(r.pe.m.nq == 0);
    REQUIRE(r.pe.m.nr == 76);
    REQUIRE(isnan(r.pe.m.sp()));
    REQUIRE(isnan(r.pe.m.sn()));

    const auto t = Test::test("-t TransAlign -rgtf data/trans/RNA.v1.gtf -ubam tests/data/trans/cufflinks.sam");

    REQUIRE(t.status == 0);
    REQUIRE(t.output.find("Transcriptome Analysis") != std::string::npos);
}

TEST_CASE("TAlign_Simulations_All_Filtered")
{
    Test::trans();
    TAlign::Options opts;
    
    for (auto i: Standard::instance().r_seqs_A)
    {
        opts.filters.insert(i.first);
    }

    const auto r = TAlign::analyze("tests/data/trans/A1/accepted_hits.sam", opts);

    REQUIRE(r.pe.m.nq == 25);
    REQUIRE(r.pe.m.nr == 76);
    REQUIRE(isnan(r.pe.m.sn()));
    REQUIRE(r.pe.m.sp() == 0.0);
    REQUIRE(r.pe.s.id == "");
    
    REQUIRE(r.pi.m.nq == 112);
    REQUIRE(r.pi.m.nr == 76);
    REQUIRE(isnan(r.pi.m.sn()));
    REQUIRE(r.pi.m.sp() == 0.0);
    REQUIRE(r.pi.s.id == "");
}

TEST_CASE("TAlign_Simulations")
{
    Test::trans();
    const auto r = TAlign::analyze("tests/data/trans/A1/accepted_hits.sam");

    REQUIRE(r.pb.m.nq == 23712);
    REQUIRE(r.pb.m.nr == 149219);
    REQUIRE(r.pb.m.sp() == Approx(0.9986082996));
    REQUIRE(r.pb.m.sn() == Approx(0.1586862263));
    REQUIRE(r.pb.s.id == "R2_24");
    REQUIRE(r.pb.s.counts == 1);
    REQUIRE(r.pb.s.abund == 1831.0);

    REQUIRE(r.pe.m.nq == 6107);
    REQUIRE(r.pe.m.nr == 6137);
    REQUIRE(r.pe.m.sp() == Approx(0.995906337));
    REQUIRE(r.pe.m.sn() == Approx(0.9910379664));
    REQUIRE(r.pe.s.id == "R2_24");
    REQUIRE(r.pe.s.counts == 1);
    REQUIRE(r.pe.s.abund == Approx(1831.0));
    
    REQUIRE(r.pi.m.nq == 1097);
    REQUIRE(r.pi.m.nr == 1047);
    REQUIRE(r.pi.m.sp() == Approx(0.8979033728));
    REQUIRE(r.pi.m.sn() == Approx(0.9407831901));
    REQUIRE(r.pi.s.id == "R2_7");
    REQUIRE(r.pi.s.counts == 1);
    REQUIRE(r.pi.s.abund == Approx(117187.0));
}