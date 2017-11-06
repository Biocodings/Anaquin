#include <catch.hpp>
#include "Kallisto.hpp"
#include <iostream>

using namespace Anaquin;

TEST_CASE("Kallisto_Shared")
{
    const auto x = KBuildIndex("tests/data/A.V.23.fa", 31);
    REQUIRE(!x.empty());

    const auto r1 = KQuery(x, "ACAGAAATAAACGAAACAGTTCTAGTAAAAA");
    const auto r2 = KQuery(x, "ATTTCCTTACTCCTCTTAAATCACCGGTGAG");
    const auto r3 = KQuery(x, "CCGGTGAGTCCGTGAACGAGTCGAGGGTAAA");
    const auto r4 = KQuery(x, "CCTTTCCCCCGTTTCCCCGACCCCTCGGAAG");
    
    REQUIRE(r1.size() == 2);
    REQUIRE(r1.count("CI_012_R"));
    REQUIRE(r1.count("CI_012_V"));
    REQUIRE(r2.size() == 2);
    REQUIRE(r2.count("CI_015_R"));
    REQUIRE(r2.count("CI_015_V"));
    REQUIRE(r3.size() == 2);
    REQUIRE(r3.count("CI_015_R"));
    REQUIRE(r3.count("CI_015_V"));
    REQUIRE(r4.size() == 2);
    REQUIRE(r4.count("CI_025_R"));
    REQUIRE(r4.count("CI_025_V"));
}

TEST_CASE("Kallisto_Unique")
{
    const auto x = KBuildIndex("tests/data/A.V.23.fa", 31);
    REQUIRE(!x.empty());
    
    const auto r1 = KQuery(x, "CAGAAATAAACGAAACAGTTCTAGTAAAAAA");
    const auto r2 = KQuery(x, "AGAAATAAACGAAACAGTTCTAGTAAAAAAC");
    const auto r3 = KQuery(x, "AATAAACGAAACAGTTCTAGTAAAAAACAAT");

    REQUIRE(r1.size() == 1);
    REQUIRE(r1.count("CI_012_R"));
    REQUIRE(r2.size() == 1);
    REQUIRE(r2.count("CI_012_R"));
    REQUIRE(r3.size() == 1);
    REQUIRE(r3.count("CI_012_R"));
}

TEST_CASE("Kallisto_2")
{
    const auto x = KBuildIndex("tests/data/A.V.23.fa", 31);
    REQUIRE(!x.empty());
    
    const auto r1 = KQuery(x, "CCTTCCCCGTCTAAAGCCCCAGATCCGAACC");
    const auto r2 = KQuery(x, "GGTTCGGATCTGGGGCTTTAGACGGGGAAGG");
    const auto r3 = KQuery(x, "GGTTCGGATCTGGGGCTTTAGACGGGGAAGA");
    const auto r4 = KQuery(x, "CCTTCCCCGTCTAAAGCCCCAGATCCGAACA");

    REQUIRE(r1.size() == 1);
    REQUIRE(r1.count("GI_030_V"));

    REQUIRE(r2.size() == 1);
    REQUIRE(r2.count("GI_030_V"));

    REQUIRE(r3.empty());
    REQUIRE(r4.empty());
}

TEST_CASE("Kallisto_3")
{
    std::map<SequinID, Base> m;
    const auto f = KHumanFA("tests/data/A.V.23.fa", m);
    const auto x = KBuildIndex(f, 31);

    const auto r1 = KQuery(x, "GTAGATCACCTGAGGTCGGGAGTTCAAGACC");
    const auto r2 = KQuery(x, "GGAAGGGGCAGATTTCGGGGTCTAGGCTTGG");
    const auto r3 = KQuery(x, "GGTTCGGATCTGGGGCTTTAGACGGGGAAGG");
    
    REQUIRE(r1.size() == 2);
    REQUIRE(r1.count("GI_044_R"));
    REQUIRE(r1.count("GI_044_V"));

    REQUIRE(r2.size() == 1);
    REQUIRE(r2.count("GI_030_V"));
    
    REQUIRE(r3.size() == 0);
}
