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

TEST_CASE("Interval_Test_3")
{
    Interval i("Test", Locus(0, 9));
    
    i.add(Locus(4, 6));
    i.add(Locus(4, 6));
    i.add(Locus(4, 6));
    i.add(Locus(1, 2));
    i.add(Locus(3, 3));
    
    std::vector<Base> x, y, z;
    
    i.bedGraph([&](const ChrID &id, Base i, Base j, Base depth)
    {
        x.push_back(i);
        y.push_back(j);
        z.push_back(depth);
    });

    REQUIRE(x.size() == 4);
    REQUIRE(y.size() == 4);
    REQUIRE(z.size() == 4);
    
    REQUIRE(x[0] == 0);
    REQUIRE(x[1] == 1);
    REQUIRE(x[2] == 4);
    REQUIRE(x[3] == 7);

    REQUIRE(y[0] == 1);
    REQUIRE(y[1] == 4);
    REQUIRE(y[2] == 7);
    REQUIRE(y[3] == 10);

    REQUIRE(z[0] == 0);
    REQUIRE(z[1] == 1);
    REQUIRE(z[2] == 3);
    REQUIRE(z[3] == 0);
}

TEST_CASE("Interval_Test_4")
{
    Intervals<> i;

    i.add(Interval("1", Locus(0,  24)));
    i.add(Interval("2", Locus(25, 49)));
    i.add(Interval("3", Locus(50, 74)));
    i.add(Interval("4", Locus(75, 99)));
    i.build();
    
    REQUIRE(i.contains(Locus(10, 20)));
    REQUIRE(i.contains(Locus(25, 26)));
    REQUIRE(i.contains(Locus(28, 28)));
    REQUIRE(i.contains(Locus(90, 95)));

    REQUIRE(!i.contains(Locus(00, 99)));
    REQUIRE(!i.contains(Locus(25, 50)));
    REQUIRE(!i.contains(Locus(51, 99)));
    REQUIRE(!i.contains(Locus(10, 60)));

    REQUIRE(i.overlap(Locus(00, 00)));
    REQUIRE(i.overlap(Locus(00, 99)));
    REQUIRE(i.overlap(Locus(25, 50)));
    REQUIRE(i.overlap(Locus(51, 99)));
    REQUIRE(i.overlap(Locus(10, 60)));
    REQUIRE(i.overlap(Locus(99, 99)));
    REQUIRE(i.overlap(Locus(0, 1000)));

    REQUIRE(!i.overlap(Locus(100,  100)));
    REQUIRE(!i.overlap(Locus(500,  600)));
    REQUIRE(!i.overlap(Locus(400,  450)));
    REQUIRE(!i.overlap(Locus(1000, 1000)));
}
