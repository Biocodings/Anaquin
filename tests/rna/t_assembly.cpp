#include <catch.hpp>
#include "rna/r_assembly.hpp"

using namespace Spike;

TEST_CASE("RAssembly_Simulation")
{
    const auto r = RAssembly::analyze("tests/data/rna/transcripts.gtf");

    REQUIRE(r.me.nq == 490);
    REQUIRE(r.me.nr == 491);
    REQUIRE(r.me.sp() == Approx(1.0));
    REQUIRE(r.me.sn() == Approx(0.9979633401));
    REQUIRE(r.se.id == "R_1_4_R");
    REQUIRE(r.se.counts == 1);
    REQUIRE(r.se.abund == Approx(0.3051757812));
    
    REQUIRE(r.mt.nq == 61);
    REQUIRE(r.mt.nr == 61);
    REQUIRE(r.mt.sp() == Approx(1.0));
    REQUIRE(r.mt.sn() == Approx(1.0));
    REQUIRE(r.st.id == "R_9_2_R");
    REQUIRE(r.st.counts == 1);
    REQUIRE(r.st.abund == Approx(0.0190734863));
    
    REQUIRE(r.mi.nq == 429);
    REQUIRE(r.mi.nr == 431);
    REQUIRE(r.mi.sp() == Approx(1.0));
    REQUIRE(r.mi.sn() == Approx(0.9953596288));
    REQUIRE(r.si.id == "R_7_2_R");
    REQUIRE(r.si.counts == 1);
    REQUIRE(r.si.abund == Approx(0.6103515625));
    
    REQUIRE(r.mb.nq == 54077);
    REQUIRE(r.mb.nr == 54280);
    REQUIRE(r.mb.sp() == Approx(1.0));
    REQUIRE(r.mb.sn() == Approx(0.9962601326));
    REQUIRE(r.sb.id == "R_1_4");
    REQUIRE(r.sb.counts == 1);
    REQUIRE(r.sb.abund == Approx(0.9155273438));
}

TEST_CASE("RAssembly_Simulations_All_Filtered")
{
    RAssembly::Options opts;
    const auto &s = Standard::instance();
    
    for (auto i: s.r_seqs_iA)
    {
        opts.filters.insert(i.first);
    }
    
    const auto r = RAssembly::analyze("tests/data/rna/transcripts.gtf", opts);
    
    REQUIRE(r.me.nq == 0);
    REQUIRE(r.me.nr == 360);
    REQUIRE(isnan(r.me.sp()));
    REQUIRE(isnan(r.me.sn()));
    REQUIRE(r.se.id == "");
    
    REQUIRE(r.mt.nq == 0);
    REQUIRE(r.mt.nr == 61);
    REQUIRE(isnan(r.mt.sp()));
    REQUIRE(isnan(r.mt.sn()));
    REQUIRE(r.st.id == "");
    
    REQUIRE(r.mi.nq == 0);
    REQUIRE(r.mi.nr == 332);
    REQUIRE(isnan(r.mi.sp()));
    REQUIRE(isnan(r.mi.sn()));
    REQUIRE(r.si.id == "");
    
    REQUIRE(r.mb.nq == 0);
    REQUIRE(r.mb.nr == 54280);
    REQUIRE(isnan(r.mb.sp()));
    REQUIRE(isnan(r.mb.sn()));
}

TEST_CASE("RAssembly_Simulations_Denovo")
{
    //const auto r = RAssembly::analyze("tests/data/rna/transcripts_dn.gtf");

    //REQUIRE(r.n  == 1040);
    //REQUIRE(r.nr == 1040);
    //REQUIRE(r.nq == 0);
}
