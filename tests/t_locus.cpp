#include "locus.hpp"
#include <catch.hpp>

using namespace Spike;

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