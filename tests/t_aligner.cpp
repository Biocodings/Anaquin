#include "gtest/gtest.h"
#include "aligner.hpp"

TEST(D1_Base_1000, AlignerTest)
{
    /*
     * Since the SAM file comes from a simulation of the reference chromosome, it's not surprising that
     * the sensitivity is 100% (from the first 1000 reads).
     */
    
    const auto stats = Aligner::base("/Users/tedwong/Sources/QA/tests/data/d1/accepted_hits.sam", Sequins(), 1000);

    ASSERT_EQ(1, stats.m.sp());
    ASSERT_TRUE(isnan(stats.m.sn()));
    ASSERT_EQ(1000, stats.m.tp);
    ASSERT_EQ(0, stats.m.fp);
    ASSERT_EQ(0, stats.m.fn);
    ASSERT_EQ(0, stats.m.tn);
    ASSERT_EQ(1000, stats.nr);
    ASSERT_EQ(1, stats.dilution);
    ASSERT_EQ(0, stats.nq);
}

TEST(D1_Splice_1000, AlignerTest)
{
    /*
     * Since the SAM file comes from a simulation of the reference chromosome, it's not surprising that
     * the sensitivity is 100% (from the first 1000 reads).
     */
    
    const auto stats = Aligner::spliced("/Users/tedwong/Sources/QA/tests/data/d1/accepted_hits.sam", Sequins(), 1000);
    
    ASSERT_EQ(1, stats.m.sp());
    ASSERT_EQ(0, stats.m.sn());
}

TEST(D2_Base, AlignerTest)
{
	// The sample file was taken from Cufflink's source distribution. It's obviously independent to the standards.
	const auto stats = Aligner::base("/Users/tedwong/Sources/QA/Tests/Data/CufflinksTest.sam");

	/*
	 * There shouldn't be any match. Both sensivity and specificity are NAN because the experiement gives no
     * power to detect anything.
	 */

    ASSERT_TRUE(isnan(stats.m.sp()));
    ASSERT_EQ(1, stats.m.sn());
    ASSERT_EQ(0, stats.m.tp);
    ASSERT_EQ(0, stats.m.fp);
    ASSERT_EQ(0, stats.m.fn);
    ASSERT_EQ(3271, stats.m.tn);
    ASSERT_EQ(0, stats.nr);
    ASSERT_EQ(0, stats.dilution);
    ASSERT_EQ(3307, stats.nq);
}