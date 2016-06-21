#include <catch.hpp>
#include "data/merged.hpp"

using namespace Anaquin;

TEST_CASE("Merged_1")
{
    MergedInterval i("Test", Locus(1, 40));
    
    Base l = 0, r = 0;
    
    i.map(Locus(3, 5), &l, &r);
    
    REQUIRE(l == 0);
    REQUIRE(r == 0);
    REQUIRE(i.size() == 1);
    
    i.map(Locus(4, 5), &l, &r);
    
    REQUIRE(l == 0);
    REQUIRE(r == 0);
    REQUIRE(i.size() == 1);
    
    i.map(Locus(4, 6), &l, &r);
    
    REQUIRE(l == 0);
    REQUIRE(r == 1);
    REQUIRE(i.size() == 1);
    
    i.map(Locus(3, 6), &l, &r);
    
    REQUIRE(l == 0);
    REQUIRE(r == 0);
    REQUIRE(i.size() == 1);
    
    i.map(Locus(8, 9), &l, &r);
    
    REQUIRE(l == 0);
    REQUIRE(r == 0);
    REQUIRE(i.size() == 2);
    
    i.map(Locus(8, 9), &l, &r);
    
    REQUIRE(l == 0);
    REQUIRE(r == 0);
    REQUIRE(i.size() == 2);
    
    i.map(Locus(1, 20), &l, &r);
    
    REQUIRE(l == 2);
    REQUIRE(r == 14); // 15?
    REQUIRE(i.size() == 1);

    i.map(Locus(1, 20), &l, &r);
    
    REQUIRE(l == 0);
    REQUIRE(r == 0);
    REQUIRE(i.size() == 1);

    const auto r1 = i.stats();
    
    REQUIRE(r1.length   == 40);
    REQUIRE(r1.nonZeros == 20);
}

TEST_CASE("Merged_2")
{
    /*
     *  (0,2) -> 0
     *  (3,5) -> 4
     *  (6,6) -> 0
     *  (7,7) -> 2
     *  (8,8) -> 1
     *  (9,9) -> 0
     *
     *  0,0,0,0,0,1,2,4,4,4
     */
    
    MergedInterval i("Test", Locus(0, 9));

    Base l = 0, r = 0;
    
    i.map(Locus(3, 5), &l, &r);

    REQUIRE(l == 0);
    REQUIRE(r == 0);
    
    i.map(Locus(3, 5), &l, &r);

    REQUIRE(l == 0);
    REQUIRE(r == 0);

    i.map(Locus(3, 5), &l, &r);

    REQUIRE(l == 0);
    REQUIRE(r == 0);

    i.map(Locus(3, 5), &l, &r);
    
    REQUIRE(l == 0);
    REQUIRE(r == 0);
    
    i.map(Locus(7, 7), &l, &r);
    
    REQUIRE(l == 0);
    REQUIRE(r == 0);
    
    i.map(Locus(7, 7), &l, &r);
    
    REQUIRE(l == 0);
    REQUIRE(r == 0);
    
    i.map(Locus(8, 8), &l, &r);

    REQUIRE(l == 0);
    REQUIRE(r == 0);
    
    const auto r1 = i.stats();
    
    REQUIRE(r1.length   == 10);
    REQUIRE(r1.nonZeros == 5);
}

TEST_CASE("Merged_3")
{
    MergedInterval i("Test", Locus(1, 40));
    
    const auto r = i.stats();
    
    REQUIRE(r.length   == 40);
    REQUIRE(r.nonZeros == 0);
}

TEST_CASE("Merged_4")
{
    MergedInterval i("Test", Locus(1, 40));
    
    i.map(Locus(1, 40));
    i.map(Locus(1, 40));
    i.map(Locus(1, 40));
    i.map(Locus(1, 40));
    i.map(Locus(1, 40));

    const auto r = i.stats();
    
    REQUIRE(r.length   == 40);
    REQUIRE(r.nonZeros == 40);
}

TEST_CASE("Merged_5")
{
    MergedInterval i("Test", Locus(1, 40));
    
    i.map(Locus(1, 20));
    i.map(Locus(1, 10));
    i.map(Locus(1, 14));
    i.map(Locus(1, 16));
    i.map(Locus(1, 19));
    
    const auto r = i.stats();
    
    REQUIRE(r.length   == 40);
    REQUIRE(r.nonZeros == 20);
}