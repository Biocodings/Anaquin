#include <catch.hpp>
#include "rna/r_align.hpp"

using namespace Spike;

static std::string exts[] = { "sam", "bam" };

TEST_CASE("RAlign_Cufflinks")
{
    // The sample file was taken from Cufflink's source distribution
    const auto r = RAlign::analyze("tests/data/cufflinks.sam");
    
    REQUIRE(0 == r.pe.m.nq);
    REQUIRE(360 == r.pe.m.nr);
    REQUIRE(isnan(r.pe.m.sp()));
    REQUIRE(isnan(r.pe.m.sn()));
}

TEST_CASE("RAlign_Simulations")
{
    for (auto ex : exts)
    {
        const auto r = RAlign::analyze("tests/data/rna/accepted_hits." + ex);

        REQUIRE(r.pb.m.nq == 31718);
        REQUIRE(r.pb.m.nr == 54280);
        REQUIRE(r.pb.m.sp() == Approx(0.9978245791));
        REQUIRE(r.pb.m.sn() == Approx(0.5830692704));
        REQUIRE(r.pb.s.id == "R_5_1");
        REQUIRE(r.pb.s.counts == 2);
        REQUIRE(r.pb.s.abund == 19.53125);
        
        REQUIRE(r.pe.m.nq == 161295);
        REQUIRE(r.pe.m.nr == 161142);
        REQUIRE(r.pe.m.sp() == Approx(0.9981400539));
        REQUIRE(r.pe.m.sn() == Approx(0.9990691613));
        REQUIRE(r.pe.s.id == "R_5_1");
        REQUIRE(r.pe.s.counts == 22);
        REQUIRE(r.pe.s.abund == Approx(19.53125));

        REQUIRE(r.pi.m.nq == 59243);
        REQUIRE(r.pi.m.nr == 59121);
        REQUIRE(r.pi.m.sp() == Approx(0.9957800922));
        REQUIRE(r.pi.m.sn() == Approx(0.9978349487));
        REQUIRE(r.pi.s.id == "R_5_1");
        REQUIRE(r.pi.s.counts == 3);
        REQUIRE(r.pi.s.abund == Approx(19.53125));
    }
}

TEST_CASE("RAlign_Simulations_All_Filtered")
{
    for (auto ex : exts)
    {
        RAlign::Options opts;
        
        for (auto i: Standard::instance().r_seqs_A)
        {
            opts.filters.insert(i.first);
        }
        
        opts.filters.insert("R_5_3_V");
        const auto r = RAlign::analyze("tests/data/rna/accepted_hits." + ex, opts);
        
        REQUIRE(r.pe.m.nq == 300);
        REQUIRE(r.pe.m.nr == 360);
        REQUIRE(isnan(r.pe.m.sn()));
        REQUIRE(r.pe.m.sp() == 0.0);
        REQUIRE(r.pe.s.id == "");

        REQUIRE(r.pi.m.nq == 250);
        REQUIRE(r.pi.m.nr == 332);
        REQUIRE(isnan(r.pi.m.sn()));
        REQUIRE(r.pi.m.sp() == 0.0);
        REQUIRE(r.pi.s.id == "");
    }
}