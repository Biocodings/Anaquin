#include <catch.hpp>
#include "rna/r_assembly.hpp"

using namespace Spike;

TEST_CASE("RAssembly_Simulation")
{
    const auto r = RAssembly::analyze("tests/data/rna/A1/transcripts.gtf");

    REQUIRE(r.pe.m.nq == 490);
    REQUIRE(r.pe.m.nr == 491);
    REQUIRE(r.pe.m.sp() == Approx(1.0));
    REQUIRE(r.pe.m.sn() == Approx(0.9979633401));
    REQUIRE(r.pe.s.id == "R_1_4_R");
    REQUIRE(r.pe.s.counts == 1);
    REQUIRE(r.pe.s.abund == Approx(0.3051757812));

    REQUIRE(r.pt.m.nq == 61);
    REQUIRE(r.pt.m.nr == 61);
    REQUIRE(r.pt.m.sp() == Approx(1.0));
    REQUIRE(r.pt.m.sn() == Approx(1.0));
    REQUIRE(r.pt.s.id == "R_9_2_R");
    REQUIRE(r.pt.s.counts == 1);
    REQUIRE(r.pt.s.abund == Approx(0.0190734863));

    REQUIRE(r.pi.m.nq == 429);
    REQUIRE(r.pi.m.nr == 431);
    REQUIRE(r.pi.m.sp() == Approx(1.0));
    REQUIRE(r.pi.m.sn() == Approx(0.9953596288));
    REQUIRE(r.pi.s.id == "R_7_2_R");
    REQUIRE(r.pi.s.counts == 1);
    REQUIRE(r.pi.s.abund == Approx(0.6103515625));
    
    REQUIRE(r.pb.m.nq == 54077);
    REQUIRE(r.pb.m.nr == 54280);
    REQUIRE(r.pb.m.sp() == Approx(1.0));
    REQUIRE(r.pb.m.sn() == Approx(0.9962601326));
    REQUIRE(r.pb.s.id == "R_1_4");
    REQUIRE(r.pb.s.counts == 1);
    REQUIRE(r.pb.s.abund == Approx(0.9155273438));
}

TEST_CASE("RAssembly_Simulations_All_Filtered")
{
    RAssembly::Options opts;
    const auto &s = Standard::instance();
    
    for (auto i: s.r_seqs_A)
    {
        opts.filters.insert(i.first);
    }
    
    const auto r = RAssembly::analyze("tests/data/rna/transcripts.gtf", opts);
    
    REQUIRE(r.pe.m.nq == 0);
    REQUIRE(r.pe.m.nr == 360);
    REQUIRE(isnan(r.pe.m.sp()));
    REQUIRE(isnan(r.pe.m.sn()));
    REQUIRE(r.pe.s.id == "");
    
    REQUIRE(r.pt.m.nq == 0);
    REQUIRE(r.pt.m.nr == 61);
    REQUIRE(isnan(r.pt.m.sp()));
    REQUIRE(isnan(r.pt.m.sn()));
    REQUIRE(r.pt.s.id == "");
    
    REQUIRE(r.pi.m.nq == 0);
    REQUIRE(r.pi.m.nr == 332);
    REQUIRE(isnan(r.pi.m.sp()));
    REQUIRE(isnan(r.pi.m.sn()));
    REQUIRE(r.pi.s.id == "");
    
    REQUIRE(r.pb.m.nq == 0);
    REQUIRE(r.pb.m.nr == 54280);
    REQUIRE(isnan(r.pb.m.sp()));
    REQUIRE(isnan(r.pb.m.sn()));
}