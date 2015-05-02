#include <catch.hpp>
#include "rna/r_align.hpp"

using namespace Spike;

static std::string exts[] = { "sam", "bam" };

TEST_CASE("RAlign_Simulations")
{
    for (auto ex : exts)
    {
        const auto r = RAlign::analyze("tests/data/rna_sims/accepted_hits." + ex);

        std::cout << r.me.nq << " " << r.me.nr << std::endl;
        
        
//        REQUIRE(155816 == r.me.n());
        REQUIRE(154637 == r.me.nq);
        REQUIRE(154638 == r.me.nr);
        
        REQUIRE(r.me.sp() == Approx(0.992382));
        REQUIRE(r.me.sn() == Approx(0.992376));
        
        REQUIRE(r.se.id == "R_5_2");
        REQUIRE(r.se.counts == 23);
        REQUIRE(r.se.abund == 9765.0);
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