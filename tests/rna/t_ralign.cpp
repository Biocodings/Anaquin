#include <catch.hpp>
#include "rna/r_align.hpp"

using namespace Spike;

static std::string exts[] = { "sam", "bam" };

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
        const auto r = RAlign::analyze("tests/data/rna_sims/accepted_hits." + ex, opts);
        
        REQUIRE(r.me.nq == 300);
        REQUIRE(r.me.nr == 161297);
        REQUIRE(isnan(r.me.sn()));
        REQUIRE(r.me.sp() == 0.0);
        REQUIRE(r.se.id == "");

        REQUIRE(r.mi.nq == 0);
        REQUIRE(r.mi.nr == 436);
        REQUIRE(isnan(r.mi.sn()));
        REQUIRE(r.mi.sp() == 0.0);
        REQUIRE(r.si.id == "");
    }
}

TEST_CASE("RAlign_Simulations")
{
    for (auto ex : exts)
    {
        const auto r = RAlign::analyze("tests/data/rna_sims/accepted_hits." + ex);

        REQUIRE(r.me.nq == 161295);
        REQUIRE(r.me.nr == 161297);
        REQUIRE(r.me.sp() == Approx(0.9981400539));
        REQUIRE(r.me.sn() == Approx(0.9981276775));
        REQUIRE(r.se.id == "R_5_1");
        REQUIRE(r.se.counts == 22);
        REQUIRE(r.se.abund == 9765.0);

        REQUIRE(r.mi.nq == 161295);
        REQUIRE(r.mi.nr == 161297);
        REQUIRE(r.mi.sp() == Approx(0.9981400539));
        REQUIRE(r.mi.sn() == Approx(0.9981276775));
        REQUIRE(r.si.id == "R_5_1");
        REQUIRE(r.si.counts == 22);
        REQUIRE(r.si.abund == 9765.0);
    }
}

TEST_CASE("RAlign_Cufflinks")
{
    // The sample file was taken from Cufflink's source distribution. It's obviously independent.
    const auto r = RAlign::analyze("tests/data/cufflinks.sam");

    REQUIRE(3329 == r.me.n());
    REQUIRE(0 == r.me.nq);
    REQUIRE(3329 == r.me.nr);
    REQUIRE(isnan(r.me.sp()));
    REQUIRE(isnan(r.me.sn()));
}