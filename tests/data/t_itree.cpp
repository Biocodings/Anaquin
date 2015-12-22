#include <catch.hpp>
#include "data/locus.hpp"
#include "data/itree.hpp"

using namespace Anaquin;

Interval_<Locus> locusToInterval(const Locus &l)
{
    return Interval_<Locus>(l.start, l.end, l);
}

TEST_CASE("ITree_Empty")
{
    IntervalTree<int> t;
    REQUIRE(t.findOverlapping(-1,1).size() == 0);
}

TEST_CASE("ITree_2")
{
    auto loci = std::vector<Interval_<Locus>>
    {
        locusToInterval(Locus(0,  24)),
        locusToInterval(Locus(25, 49)),
        locusToInterval(Locus(50, 74)),
        locusToInterval(Locus(75, 99)),
    };
    
    IntervalTree<Locus> t { loci };
    
    SECTION ("findContains")
    {
        auto v = t.findContains(10, 20);
        REQUIRE(v.size() == 1);
    }
}

TEST_CASE("ITree_1")
{
    auto loci = std::vector<Interval_<Locus>>
    {
        locusToInterval(Locus(0,  10)),
        locusToInterval(Locus(11, 20)),
        locusToInterval(Locus(21, 30)),
        locusToInterval(Locus(31, 40)),
        locusToInterval(Locus(41, 50)),
        locusToInterval(Locus(51, 60)),
        locusToInterval(Locus(61, 70)),
        locusToInterval(Locus(71, 80)),
        locusToInterval(Locus(81, 90)),
        locusToInterval(Locus(91, 100)),
    };
    
    IntervalTree<Locus> t { loci };

    SECTION ("findOverlapping")
    {
        auto v = t.findOverlapping(85, 99);
        REQUIRE(v.size() == 2);
        REQUIRE(v.front().start == 81);
        REQUIRE(v.front().stop  == 90);
        REQUIRE(v.back().start  == 91);
        REQUIRE(v.back().stop   == 100);
    }

    SECTION ("findContained")
    {
        auto v = t.findContained(85, 99);
        REQUIRE(v.size() == 0);
    }
}