#include <vector>
#include "gtest/gtest.h"
#include "StandardFactory.hpp"

using namespace std;

TEST(TestID, StandardFactoryTest)
{
    const auto r = StandardFactory::reference();
	ASSERT_EQ("chrT", r.id);
}

TEST(TestJunction, StandardFactoryTest)
{
    const auto r  = StandardFactory::reference();
    const auto &g = r.genes[0];
    
    ASSERT_EQ(9, g.exons.size());
    ASSERT_EQ(8, g.js.size());

    ASSERT_EQ(388626, g.js[0].loc.start());
    ASSERT_EQ(447005, g.js[0].loc.end());

    ASSERT_EQ(447068, g.js[1].loc.start());
    ASSERT_EQ(480183, g.js[1].loc.end());

    ASSERT_EQ(480429, g.js[2].loc.start());
    ASSERT_EQ(524943, g.js[2].loc.end());
    
    ASSERT_EQ(525182, g.js[3].loc.start());
    ASSERT_EQ(527936, g.js[3].loc.end());

    ASSERT_EQ(528162, g.js[4].loc.start());
    ASSERT_EQ(538413, g.js[4].loc.end());

    ASSERT_EQ(538493, g.js[5].loc.start());
    ASSERT_EQ(539348, g.js[5].loc.end());

    ASSERT_EQ(539742, g.js[6].loc.start());
    ASSERT_EQ(540365, g.js[6].loc.end());

    ASSERT_EQ(540493, g.js[7].loc.start());
    ASSERT_EQ(542286, g.js[7].loc.end());
}

TEST(TestSequins, StandardFactoryTest)
{
	const auto r = StandardFactory::reference();
	ASSERT_EQ(62, r.sequins.size());

    const auto ids = { "R_1_1_R",  "R_1_1_V",  "R_1_2_R",
                       "R_1_2_V",  "R_1_3_R",  "R_1_3_V",
                       "R_1_4_R",  "R_10_1_R", "R_10_1_V",
                       "R_10_2_R", "R_10_2_V", "R_10_3_R",
                       "R_10_3_V", "R_2_1_R",  "R_2_1_V",
                       "R_2_2_R",  "R_2_2_V",  "R_2_3_R",
                       "R_2_3_V",  "R_2_4_R",  "R_3_1_R",
                       "R_3_1_V",  "R_3_2_R",  "R_3_2_V",
                       "R_3_3_R",  "R_3_3_V",  "R_4_1_R",
                       "R_4_1_V",  "R_4_2_R",  "R_4_2_V",
                       "R_4_3_R",  "R_4_3_V",  "R_5_1_R",
                       "R_5_1_V",  "R_5_2_R",  "R_5_2_V",
                       "R_5_3_R",  "R_5_3_V",  "R_6_1_R",
                       "R_6_1_V",  "R_6_2_R",  "R_6_2_V",
                       "R_6_3_R",  "R_6_3_V",  "R_7_1_R",
                       "R_7_1_V",  "R_7_2_R",  "R_7_2_V",
                       "R_7_3_R",  "R_7_3_V",  "R_8_1_R",
                       "R_8_1_V",  "R_8_2_R",  "R_8_2_V",
                       "R_8_3_R",  "R_8_3_V",  "R_9_1_R",
                       "R_9_1_V",  "R_9_2_R",  "R_9_2_V",
                       "R_9_3_R",  "R_9_3_V"
                    };

    for (auto id : ids)
    {
        ASSERT_TRUE(r.sequins.count(id));
    }    
}

