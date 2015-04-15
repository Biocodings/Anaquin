#include <vector>
#include <catch.hpp>
#include "standard_factory.hpp"

using namespace Spike;

TEST_CASE("Standard_ID")
{
    const auto r = StandardFactory::reference();
	REQUIRE("chrT" == r.id);
}

TEST_CASE("Standard_Isoforms")
{
    const auto r = StandardFactory::reference();
    
    REQUIRE(r.seqs_iA.size() == 62);
    REQUIRE(r.seqs_iB.size() == 62);

    const auto ids =
    {
        "R_1_1_R",  "R_1_1_V",  "R_1_2_R",  "R_1_2_V",  "R_1_3_R",  "R_1_3_V",
        "R_1_4_R",  "R_10_1_R", "R_10_1_V", "R_10_2_R", "R_10_2_V", "R_10_3_R",
        "R_10_3_V", "R_2_1_R",  "R_2_1_V",  "R_2_2_R",  "R_2_2_V",  "R_2_3_R",
        "R_2_3_V",  "R_2_4_R",  "R_3_1_R",  "R_3_1_V",  "R_3_2_R",  "R_3_2_V",
        "R_3_3_R",  "R_3_3_V",  "R_4_1_R",  "R_4_1_V",  "R_4_2_R",  "R_4_2_V",
        "R_4_3_R",  "R_4_3_V",  "R_5_1_R",  "R_5_1_V",  "R_5_2_R",  "R_5_2_V",
        "R_5_3_R",  "R_6_1_R",  "R_6_1_V",  "R_6_2_R",  "R_6_2_V",  "R_6_3_R",
        "R_6_3_V",  "R_7_1_R",  "R_7_1_V",  "R_7_2_R",  "R_7_2_V",  "R_7_3_R",
        "R_7_3_V",  "R_8_1_R",  "R_8_1_V",  "R_8_2_R",  "R_8_2_V",  "R_8_3_R",
        "R_8_3_V",  "R_9_1_R",  "R_9_1_V",  "R_9_2_R",  "R_9_2_V",  "R_9_3_R",
        "R_9_3_V"
    };

    for (auto id : ids)
    {
        REQUIRE(r.seqs_iA.count(id));
    }
}