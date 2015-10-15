#include <catch.hpp>
#include "data/intervals.hpp"

using namespace Anaquin;

TEST_CASE("Interval_Test_1")
{
    /*
     * (0,2) -> 0
     * (3,5) -> 4
     * (6,6) -> 0
     * (7,7) -> 2
     * (8,8) -> 1
     * (9,9) -> 0
     *
     *  0,0,0,0,0,1,2,4,4,4
     */
    
    Interval i("Test", Locus(0, 9));

    i.add(Locus(3, 5));
    i.add(Locus(3, 5));
    i.add(Locus(3, 5));
    i.add(Locus(3, 5));
    i.add(Locus(7, 7));
    i.add(Locus(7, 7));
    i.add(Locus(8, 8));

    const auto r = i.stats();
    
    REQUIRE(r.min      == 0);
    REQUIRE(r.max      == 4);
    REQUIRE(r.length   == 10);
    REQUIRE(r.zeros    == 5);
    REQUIRE(r.nonZeros == 5);
    REQUIRE(r.sums     == 15);
    REQUIRE(r.p25      == 0);
    REQUIRE(r.p50      == 0);
    REQUIRE(r.p75      == 2);
}

TEST_CASE("Interval_Test_2")
{
    /*
     * (0,4) -> 5
     * (5,5) -> 1
     * (6,9) -> 5
     *
     *  1,5,5,5,5,5,5,5,5
     */
    
    Interval i("Test", Locus(0, 9));
    
    i.add(Locus(0, 4));
    i.add(Locus(0, 4));
    i.add(Locus(0, 4));
    i.add(Locus(0, 4));
    i.add(Locus(0, 4));
    i.add(Locus(5, 5));
    i.add(Locus(6, 9));
    i.add(Locus(6, 9));
    i.add(Locus(6, 9));
    i.add(Locus(6, 9));
    i.add(Locus(6, 9));

    const auto r = i.stats();
    
    REQUIRE(r.min      == 1);
    REQUIRE(r.max      == 5);
    REQUIRE(r.length   == 10);
    REQUIRE(r.zeros    == 0);
    REQUIRE(r.nonZeros == 10);
    REQUIRE(r.sums     == 46);
    REQUIRE(r.p25      == 5);
    REQUIRE(r.p50      == 5);
    REQUIRE(r.p75      == 5);
}
    


