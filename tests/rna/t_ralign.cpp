#include <catch.hpp>
#include "rna/r_align.hpp"

using namespace Spike;

static std::string exts[] = { "sam", "bam" };

#ifdef REGRESSION_TESTING

TEST_CASE("RAlign_RK1A")
{
    const auto r = RAlign::analyze("tests/data/R_K562/RK1A_filtered.bam");
}

TEST_CASE("RAlign_RK2A")
{
    const auto r = RAlign::analyze("tests/data/R_K562/RK2A_filtered.bam");
}

TEST_CASE("RAlign_RK3A")
{
    const auto r = RAlign::analyze("tests/data/R_K562/RK3A_filtered.bam");
}

#endif

TEST_CASE("RAlign_Simulations")
{
    for (auto ex : exts)
    {
        const auto r = RAlign::analyze("tests/data/rna/accepted_hits." + ex);

        REQUIRE(r.mb.nq == 31718);
        REQUIRE(r.mb.nr == 54280);
        REQUIRE(r.mb.sp() == Approx(0.9978245791));
        REQUIRE(r.mb.sn() == Approx(0.5830692704));
        REQUIRE(r.sb.id == "R_5_1");
        REQUIRE(r.sb.counts == 2);
        REQUIRE(r.sb.abund == 9.765625);
        
        REQUIRE(r.me.nq == 161295);
        REQUIRE(r.me.nr == 161142);
        REQUIRE(r.me.sp() == Approx(0.9981400539));
        REQUIRE(r.me.sn() == Approx(0.9990691613));
        REQUIRE(r.se.id == "R_5_1");
        REQUIRE(r.se.counts == 22);
        REQUIRE(r.se.abund == Approx(9.765625));

        REQUIRE(r.mi.nq == 59243);
        REQUIRE(r.mi.nr == 59121);
        REQUIRE(r.mi.sp() == Approx(0.9957800922));
        REQUIRE(r.mi.sn() == Approx(0.9978349487));
        REQUIRE(r.si.id == "R_5_1");
        REQUIRE(r.si.counts == 3);
        REQUIRE(r.si.abund == Approx(9.765625));
    }
}

TEST_CASE("RAlign_Simulations_All_Filtered")
{
    for (auto ex : exts)
    {
        RAlign::Options opts;
        
        for (auto i: Standard::instance().r_seqs_iA)
        {
            opts.filters.insert(i.first);
        }
        
        opts.filters.insert("R_5_3_V");
        const auto r = RAlign::analyze("tests/data/rna/accepted_hits." + ex, opts);
        
        REQUIRE(r.me.nq == 300);
        REQUIRE(r.me.nr == 360);
        REQUIRE(isnan(r.me.sn()));
        REQUIRE(r.me.sp() == 0.0);
        REQUIRE(r.se.id == "");

        REQUIRE(r.mi.nq == 250);
        REQUIRE(r.mi.nr == 332);
        REQUIRE(isnan(r.mi.sn()));
        REQUIRE(r.mi.sp() == 0.0);
        REQUIRE(r.si.id == "");
    }
}

TEST_CASE("RAlign_Cufflinks")
{
    // The sample file was taken from Cufflink's source distribution. It's obviously independent.
    const auto r = RAlign::analyze("tests/data/cufflinks.sam");

    REQUIRE(0 == r.me.nq);
    REQUIRE(360 == r.me.nr);
    REQUIRE(isnan(r.me.sp()));
    REQUIRE(isnan(r.me.sn()));
}