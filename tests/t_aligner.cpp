#include <catch.hpp>
#include "aligner.hpp"

TEST_CASE("D1_Base_1000")
{
    /*
     * Since the SAM file comes from a simulation of the reference chromosome, it's not surprising that
     * the sensitivity is 100% (from the first 1000 reads).
     */
    
    const auto stats = Aligner::base("tests/data/d1/accepted_hits.sam", Sequins(), 1000);

    REQUIRE(1 == stats.m.sp());
    REQUIRE(isnan(stats.m.sn()));
    REQUIRE(1000 == stats.m.tp);
    REQUIRE(0 == stats.m.fp);
    REQUIRE(0 == stats.m.fn);
    REQUIRE(0 == stats.m.tn);
    REQUIRE(1000 == stats.nr);
    REQUIRE(1 == stats.dilution);
    REQUIRE(0 == stats.nq);
}

TEST_CASE("D1_Splice_1000")
{
    /*
     * Since the SAM file comes from a simulation of the reference chromosome, it's not surprising that
     * the sensitivity is 100% (from the first 1000 reads).
     */
    
    const auto stats = Aligner::spliced("tests/data/d1/accepted_hits.sam", Sequins(), 1000);
    
    REQUIRE(1 == stats.m.sp());
    REQUIRE(0 == stats.m.sn());
}

TEST_CASE("D2_Base")
{
	// The sample file was taken from Cufflink's source distribution. It's obviously independent to the standards.
	const auto stats = Aligner::base("tests/Data/CufflinksTest.sam");

	/*
	 * There shouldn't be any match. Both sensivity and specificity are NAN because the experiement gives no
     * power to detect anything.
	 */

    REQUIRE(isnan(stats.m.sp()));
    REQUIRE(1 == stats.m.sn());
    REQUIRE(0 == stats.m.tp);
    REQUIRE(0 == stats.m.fp);
    REQUIRE(0 == stats.m.fn);
    REQUIRE(3271 == stats.m.tn);
    REQUIRE(0 == stats.nr);
    REQUIRE(0 == stats.dilution);
    REQUIRE(3307 == stats.nq);
}