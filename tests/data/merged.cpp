#include <catch.hpp>
#include "data/merged.hpp"

using namespace Anaquin;

//inline MergedIntervals<> abcd(const MergedIntervals<> &data)
//{
//    const Locus *curr;
//    const Locus *next;
//
//    std::vector<const MergedInterval *> ptrs;
//
//    for (const auto &i : data._inters)
//    {
//        ptrs.push_back(&i.second);
//    }
//
//    Locus l;
//    MergedIntervals<> r;
//
//    Base a, b, x, y;
//
//    a = ptrs[0]->l().start;
//    b = ptrs[0]->l().end;
//
//    for (auto i = 1; i < ptrs.size(); i++)
//    {
//        x = ptrs[i]->l().start;
//        y = ptrs[i]->l().end;
//
//        next = (i == ptrs.size()-1) ? nullptr : &(ptrs[i+1]->l());
//
//        // The next region is non-overlapping
//        if (x > b)
//        {
//            r.add(MergedInterval(std::to_string(i), Locus(a, b)));
//            a = x;
//            b = y;
//        }
//        
//        // The next region is overlapping but outside
//        else if (y > b)
//        {
//            r.add(MergedInterval(std::to_string(i), Locus(a, x-1)));
//            a = x;
//        }
//        
//        // The next region is inside
//        else
//        {
//            r.add(MergedInterval(std::to_string(i), Locus(a, x-1)));
//            a = x;
//        }
//    }
//    
//    return r;
//}
//
//TEST_CASE("Merged_Bins")
//{
//    /*
//     *   (1,4)
//     *   (5,10)
//     *   (10,13)
//     *   (13,15)
//     *   (16,20)
//     */
//    
//    MergedIntervals<> x;
//    
//    x.add(MergedInterval("1", Locus(1,  10)));
//    x.add(MergedInterval("2", Locus(5,  15)));
//    x.add(MergedInterval("3", Locus(13, 20)));
//    x.add(MergedInterval("4", Locus(30, 40)));
//    x.add(MergedInterval("5", Locus(50, 70)));
//    x.add(MergedInterval("6", Locus(55, 80)));
//    
//    const auto r = abcd(x);
//
//    REQUIRE(r.size() == 10);
//}

using namespace Anaquin;

TEST_CASE("Merged_11")
{
    MergedInterval i("Test", Locus(1, 1000));
    
    i.map(Locus(2, 10));

    REQUIRE(i.size() == 1);
    REQUIRE(i._data.count(10));
    REQUIRE(!i._data.count(15));
    REQUIRE(i._data[10].start == 2);
    REQUIRE(i._data[10].end == 10);
    
    i.map(Locus(5, 15));
    
    REQUIRE(i.size() == 1);
    REQUIRE(!i._data.count(10));
    REQUIRE(i._data.count(15));
    REQUIRE(i._data[15].start == 2);
    REQUIRE(i._data[15].end == 15);
}

TEST_CASE("Merged_10")
{
    MergedInterval i("Test", Locus(1, 50000));
    
    i.map(Locus(2,115));
    i.map(Locus(167,283));
    i.map(Locus(334,434));
    i.map(Locus(463,563));
    i.map(Locus(618,804));
    i.map(Locus(956,1056));
    i.map(Locus(1088,1188));
    i.map(Locus(1348,1448));
    i.map(Locus(1642,1742));
    i.map(Locus(1831,1931));
    i.map(Locus(2002,2211));
    i.map(Locus(2293,2393));
    i.map(Locus(2806,2994));
    i.map(Locus(3049,3105));
    i.map(Locus(3347,3447));
    i.map(Locus(3563,3663));
    i.map(Locus(3774,3874));
    i.map(Locus(4030,4155));
    i.map(Locus(4330,4430));
    i.map(Locus(4602,4722));
    i.map(Locus(4752,4936));
    i.map(Locus(4979,5144));
    i.map(Locus(5407,5507));
    i.map(Locus(5583,5586));
    i.map(Locus(5667,5767));
    i.map(Locus(5869,6094));
    i.map(Locus(6173,6342));
    i.map(Locus(6610,6710));
    i.map(Locus(6799,6899));
    i.map(Locus(6939,7052));
    i.map(Locus(7113,7215));
    i.map(Locus(7329,7429));
    i.map(Locus(7688,7788));
    i.map(Locus(7839,7939));
    i.map(Locus(8600,8700));
    i.map(Locus(8795,8895));
    i.map(Locus(9107,9243));
    i.map(Locus(9306,9489));
    i.map(Locus(9506,9606));
    i.map(Locus(9705,9825));
    i.map(Locus(9902,10002));
    i.map(Locus(10262,10404));
    i.map(Locus(10466,10581));
    i.map(Locus(10724,10824));
    i.map(Locus(10830,10930));
    i.map(Locus(10945,11045));
    i.map(Locus(11059,11159));
    i.map(Locus(11196,11296));
    i.map(Locus(11489,11589));
    i.map(Locus(12458,12558));
    i.map(Locus(12633,12733));
    i.map(Locus(12761,12861));
    i.map(Locus(12944,13094));
    i.map(Locus(13108,13302));
    i.map(Locus(13315,13415));
    i.map(Locus(13466,13588));
    i.map(Locus(13603,13756));
    i.map(Locus(13948,14048));
    i.map(Locus(14086,14186));
    i.map(Locus(14213,14313));
    i.map(Locus(15184,15284));
    
    /*
     * i.map(Locus(2002,2211));
     * i.map(Locus(2293,2393));
     */
    
    i.map(Locus(2180, 2280));
    //REQUIRE(i.size() == 2);
}

TEST_CASE("Merged_8")
{
    MergedInterval i("Test", Locus(1, 1000));
    
    i.map(Locus(2, 10));
    i.map(Locus(500, 600));
    
    REQUIRE(i.size() == 2);
    
    i.map(Locus(2, 20));
    REQUIRE(i.size() == 2);
}

TEST_CASE("Merged_9")
{
    MergedInterval i("Test", Locus(1, 1000));
    
    i.map(Locus(2, 10));
    i.map(Locus(500, 600));
    
    REQUIRE(i.size() == 2);
    
    i.map(Locus(100, 200));
    REQUIRE(i.size() == 3);
}

TEST_CASE("Merged_7")
{
    MergedInterval i("Test", Locus(1, 1000));
    
    i.map(Locus(2, 10));
    i.map(Locus(500, 600));
    
    REQUIRE(i.size() == 2);
    
    i.map(Locus(172, 309));
    REQUIRE(i.size() == 3);
}

TEST_CASE("Merged_6")
{
    MergedInterval i("Test", Locus(1, 1000));

    i.map(Locus(2, 283));
    i.map(Locus(500, 600));
    
    REQUIRE(i.size() == 2);
    
    i.map(Locus(172, 309));
    REQUIRE(i.size() == 2);
}

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