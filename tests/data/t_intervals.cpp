#include <catch.hpp>
#include "data/intervals.hpp"

using namespace Anaquin;

TEST_CASE("Interval_Test_5")
{
    /*
     * Simulating a genome. Genome can be represented by multiple intervals. Thus, if
     * we have four genes, we'll have four intervals. The representation allows flexibility
     * and extensibility.
     */
    
    std::map<ChrID, Intervals<>> i;
    
    i["chrIS_A"].add(Interval("Gene A", Locus(0, 999)));
    i["chrIS_B"].add(Interval("Gene B", Locus(0, 999)));
    i["chrIS_C"].add(Interval("Gene C", Locus(0, 999)));
    i["chrIS_D"].add(Interval("Gene D", Locus(0, 999)));
    
    i["chrIS_A"].build();
    i["chrIS_B"].build();
    i["chrIS_C"].build();
    i["chrIS_D"].build();
    
    REQUIRE(i["chrIS_A"].find("Gene A"));
    REQUIRE(i["chrIS_B"].find("Gene B"));
    REQUIRE(i["chrIS_C"].find("Gene C"));
    REQUIRE(i["chrIS_D"].find("Gene D"));
    
    REQUIRE(!i["chrIS_A"].find("gene a"));
    REQUIRE(!i["chrIS_B"].find("gene b"));
    REQUIRE(!i["chrIS_C"].find("gene c"));
    REQUIRE(!i["chrIS_D"].find("gene d"));
    
    REQUIRE(i["chrIS_A"].find("Gene A")->l().length() == 1000);
    REQUIRE(i["chrIS_B"].find("Gene B")->l().length() == 1000);
    REQUIRE(i["chrIS_C"].find("Gene C")->l().length() == 1000);
    REQUIRE(i["chrIS_D"].find("Gene D")->l().length() == 1000);
    
    REQUIRE(i["chrIS_A"].contains(Locus(0, 20)));
    REQUIRE(i["chrIS_B"].contains(Locus(0, 20)));
    REQUIRE(i["chrIS_C"].contains(Locus(0, 20)));
    REQUIRE(i["chrIS_D"].contains(Locus(0, 20)));
    
    REQUIRE(!i["chrIS_A"].contains(Locus(0, 1000)));
    REQUIRE(!i["chrIS_B"].contains(Locus(0, 1000)));
    REQUIRE(!i["chrIS_C"].contains(Locus(0, 1000)));
    REQUIRE(!i["chrIS_D"].contains(Locus(0, 1000)));

    REQUIRE(i["chrIS_A"].stats().covered() == 0);
    REQUIRE(i["chrIS_B"].stats().covered() == 0);
    REQUIRE(i["chrIS_C"].stats().covered() == 0);
    REQUIRE(i["chrIS_D"].stats().covered() == 0);

    i["chrIS_A"].find("Gene A")->add(Locus(5, 14));
    i["chrIS_B"].find("Gene B")->add(Locus(5, 14));
    i["chrIS_C"].find("Gene C")->add(Locus(5, 14));
    i["chrIS_D"].find("Gene D")->add(Locus(5, 14));
 
    REQUIRE(i["chrIS_A"].stats().covered() == 0.01);
    REQUIRE(i["chrIS_B"].stats().covered() == 0.01);
    REQUIRE(i["chrIS_C"].stats().covered() == 0.01);
    REQUIRE(i["chrIS_D"].stats().covered() == 0.01);

    i["chrIS_A"].find("Gene A")->add(Locus(0, 999));
    i["chrIS_B"].find("Gene B")->add(Locus(0, 999));
    i["chrIS_C"].find("Gene C")->add(Locus(0, 999));
    i["chrIS_D"].find("Gene D")->add(Locus(0, 999));

    REQUIRE(i["chrIS_A"].stats().covered() == 1.0);
    REQUIRE(i["chrIS_B"].stats().covered() == 1.0);
    REQUIRE(i["chrIS_C"].stats().covered() == 1.0);
    REQUIRE(i["chrIS_D"].stats().covered() == 1.0);
}

TEST_CASE("Interval_Test_1")
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
    
    Interval i("Test", Locus(0, 9));

    const auto r1 = i.stats();

    REQUIRE(r1.length == 10);
    
    i.add(Locus(3, 5));
    i.add(Locus(3, 5));
    i.add(Locus(3, 5));
    i.add(Locus(3, 5));
    i.add(Locus(7, 7));
    i.add(Locus(7, 7));
    i.add(Locus(8, 8));

    const auto r2 = i.stats();
    
    REQUIRE(r2.min      == 0);
    REQUIRE(r2.max      == 4);
    REQUIRE(r2.length   == 10);
    REQUIRE(r2.zeros    == 5);
    REQUIRE(r2.nonZeros == 5);
    REQUIRE(r2.sums     == 15);
    REQUIRE(r2.p25      == 0);
    REQUIRE(r2.p50      == 0);
    REQUIRE(r2.p75      == 2);
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