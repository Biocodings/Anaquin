#include "gtest/gtest.h"
#include "AlignerAnalyst.hpp"

using namespace std;

TEST(Cufflink, AlignerStatsTest)
{
	// The sample file was taken from Cufflink's source code. It's obviously independent to the standards.
	const auto stats = AlignerAnalyst::analyze("../Tests/Data/CufflinksTest.sam");

	/*
	 * There shouldn't be any match to the chromosome. Sensivitiy is NAN because there is no data
	 * to detect whethetr the experiment is able to detect anything postively. Specificity is one
	 * because none of the reads that fails to map to the chromosome comes from it.
	 */

	ASSERT_TRUE(isnan(stats.sp));
	ASSERT_EQ(1, stats.sn);
	ASSERT_EQ(0, stats.tp);
	ASSERT_EQ(0, stats.fp);
	ASSERT_EQ(0, stats.fn);
	ASSERT_EQ(36, stats.tn);
	ASSERT_EQ(0, stats.n_chromo);
	ASSERT_EQ(3307, stats.n_sample);
}