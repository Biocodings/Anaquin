#include <math.h>
#include "gtest/gtest.h"
#include "AlignerAnalyst.hpp"

using namespace std;

TEST(Cufflink, AlignerStatsTest)
{
	// The sample file was taken from Cufflink's source distribution. It's obviously independent to the standards.
	const auto stats = AlignerAnalyst::analyze("/Users/user1/Sources/QA/Tests/Data/CufflinksTest.sam");

	/*
	 * There shouldn't be any match. Both sensivity and specificity are NAN because the experiement gives no
     * power to detect anything.
	 */

	ASSERT_TRUE(isnan(stats.sp));
    ASSERT_TRUE(isnan(stats.sn));
	ASSERT_EQ(0, stats.tp);
	ASSERT_EQ(0, stats.fp);
	ASSERT_EQ(0, stats.fn);
	ASSERT_EQ(0, stats.n_r);
	ASSERT_EQ(3307, stats.n_q);
}