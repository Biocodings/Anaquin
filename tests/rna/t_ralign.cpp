#include <catch.hpp>
#include "rna/r_align.hpp"

using namespace Spike;

TEST_CASE("RAlign_Cufflinks")
{
    // The sample file was taken from Cufflink's source distribution
    const auto r = RAlign::analyze("tests/data/rna/cufflinks.sam");
    
    REQUIRE(r.pe.m.nq == 0);
    REQUIRE(r.pe.m.nr == 78);
    REQUIRE(isnan(r.pe.m.sp()));
    REQUIRE(isnan(r.pe.m.sn()));
}

TEST_CASE("RAlign_Simulations_All_Filtered")
{
    RAlign::Options opts;
    
    for (auto i: Standard::instance().r_seqs_A)
    {
        opts.filters.insert(i.first);
    }

    const auto r = RAlign::analyze("tests/data/rna/A1/accepted_hits.sam", opts);

    REQUIRE(r.pe.m.nq == 25);
    REQUIRE(r.pe.m.nr == 78);
    REQUIRE(isnan(r.pe.m.sn()));
    REQUIRE(r.pe.m.sp() == 0.0);
    REQUIRE(r.pe.s.id == "");
    
    REQUIRE(r.pi.m.nq == 112);
    REQUIRE(r.pi.m.nr == 78);
    REQUIRE(isnan(r.pi.m.sn()));
    REQUIRE(r.pi.m.sp() == 0.0);
    REQUIRE(r.pi.s.id == "");
}

TEST_CASE("RAlign_Simulations")
{
    const auto r = RAlign::analyze("tests/data/rna/A1/accepted_hits.sam");

    REQUIRE(r.pb.m.nq == 23712);
    REQUIRE(r.pb.m.nr == 151192);
    REQUIRE(r.pb.m.sp() == Approx(0.9986082996));
    REQUIRE(r.pb.m.sn() == Approx(0.1566154294));
    REQUIRE(r.pb.s.id == "R1_22");
    REQUIRE(r.pb.s.counts == 1);
    REQUIRE(r.pb.s.abund == 10.0);
    
    REQUIRE(r.pe.m.nq == 6107);
    REQUIRE(r.pe.m.nr == 6139);
    REQUIRE(r.pe.m.sp() == Approx(0.995906337));
    REQUIRE(r.pe.m.sn() == Approx(0.9907151002));
    REQUIRE(r.pe.s.id == "R1_22");
    REQUIRE(r.pe.s.counts == 1);
    REQUIRE(r.pe.s.abund == Approx(10.0));
    
    REQUIRE(r.pi.m.nq == 1097);
    REQUIRE(r.pi.m.nr == 1049);
    REQUIRE(r.pi.m.sp() == Approx(0.8979033728));
    REQUIRE(r.pi.m.sn() == Approx(0.9389895138));
    REQUIRE(r.pi.s.id == "R2_7");
    REQUIRE(r.pi.s.counts == 1);
    REQUIRE(r.pi.s.abund == Approx(10.0));
}