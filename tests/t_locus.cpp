#include "locus.hpp"
#include <catch.hpp>

using namespace Spike;

TEST_CASE("Locus_Merge_None_Overlap")
{
    const auto m = Locus::merge(std::vector<Locus> { Locus(10, 20), Locus(21, 30), Locus(31, 40) });

    REQUIRE(m.size() == 3);
    REQUIRE(m[0] == Locus(10, 20));
    REQUIRE(m[1] == Locus(21, 30));
    REQUIRE(m[2] == Locus(31, 40));
}

TEST_CASE("Locus_Merge_Some_Overlap")
{
    const auto m = Locus::merge(std::vector<Locus> { Locus(10, 20), Locus(15, 22), Locus(31, 40) });
    
    REQUIRE(m.size() == 2);
    REQUIRE(m[0] == Locus(10, 22));
    REQUIRE(m[1] == Locus(31, 40));
}

TEST_CASE("Locus_Merge_All_Overlap")
{
    const auto m = Locus::merge(std::vector<Locus> { Locus(10, 20), Locus(10, 20), Locus(10, 20) });
    
    REQUIRE(m.size() == 1);
    REQUIRE(m[0] == Locus(10, 20));
}

TEST_CASE("Locus_No_Overlap")
{
    Locus s1(1, 100);
    Locus s2(101, 200);
    
    REQUIRE(s1.overlap(s2) == 0);
    REQUIRE(s2.overlap(s1) == 0);
}

TEST_CASE("Local_Contained")
{
    Locus s1(10, 1000);
    Locus s2(50, 60);

    REQUIRE(s1.overlap(s2) == 11);
    REQUIRE(s2.overlap(s1) == 11);
}

TEST_CASE("Locus_Overalap_One")
{
    Locus s1(1, 101);
    Locus s2(101, 200);

    REQUIRE(s1.overlap(s2) == 1);
    REQUIRE(s2.overlap(s1) == 1);
}

TEST_CASE("Locus_Equal")
{
    Locus s1(1, 100);
    REQUIRE(s1.overlap(s1) == 100);
}

TEST_CASE("Locus_Overlap")
{
    Locus s1(1, 10);
    Locus s2(5, 15);

    REQUIRE(s1.overlap(s2) == 6);
    REQUIRE(s2.overlap(s1) == 6);
}