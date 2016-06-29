#include <catch.hpp>
#include "parsers/parser_sam.hpp"

using namespace Anaquin;

TEST_CASE("Test_Insert")
{
    std::map<std::size_t, std::vector<Locus>> r;
    
    ParserSAM::parse("tests/data/insert.sam", [&](ParserSAM::Data &x, const ParserSAM::Info &)
    {
        Locus l;
        bool spliced;
        auto i = r.size();
        
        while (x.nextCigar(l, spliced))
        {
            r[i].push_back(l);
        }
    });
    
    /*
     * 8288747	60	23M2I58M42S
     */
    
    REQUIRE(r.size() == 1);
    REQUIRE(r[0].size() == 2);
    
    REQUIRE(r[0][0].start == 8288747);
    REQUIRE(r[0][0].end   == 8288769);
    REQUIRE(r[0][1].start == 8288770);
    REQUIRE(r[0][1].end   == 8288827);
}

TEST_CASE("Test_Deletion")
{
    std::vector<ParserSAM::Data> r1;
    std::vector<ParserSAM::Info> r2;
    std::map<std::size_t, std::vector<Locus>> r3;

    ParserSAM::parse("tests/data/deletion.sam", [&](ParserSAM::Data &x, const ParserSAM::Info &i)
    {
        Locus l;
        bool spliced;
        
        while (x.nextCigar(l, spliced))
        {
            r3[r1.size()].push_back(l);
        }
        
        r1.push_back(x);
        r2.push_back(i);
    });
    
    REQUIRE(r1.size() == 2);
    REQUIRE(r2.size() == 2);
    
    /*
     * 7058781	60	58M9D67M
     */
    
    REQUIRE(r3[1].size() == 2);

    REQUIRE(r3[1][0].start == 7058781);
    REQUIRE(r3[1][0].end == 7058838);

    REQUIRE(r3[1][1].start == 7058848);
    REQUIRE(r3[1][1].end == 7058914);
}

TEST_CASE("Test_SoftClip")
{
    std::vector<ParserSAM::Data> r1;
    std::vector<ParserSAM::Info> r2;
    
    ParserSAM::parse("tests/data/clip.sam", [&](const ParserSAM::Data &x, const ParserSAM::Info &i)
    {
        r1.push_back(x);
        r2.push_back(i);
    });
    
    REQUIRE(r1.size() == 2);
    REQUIRE(r2.size() == 2);

    REQUIRE(r2[0].clip);
    REQUIRE(r2[1].clip);
    
    /*
     * 4106432	60	34M76S	=
     */

    REQUIRE(r1[0].l.start == 4106432);
    REQUIRE(r1[0].l.end   == 4106465);
    
    /*
     * 4106431	60	80S35M	=
     */
    
    REQUIRE(r1[1].l.start == 4106431);
    REQUIRE(r1[1].l.end   == 4106465);
}

//TEST_CASE("Test_Junction")
//{
//    std::vector<Alignment> aligns;
//    
//    /*
//     * Cigar: 388598 - 29M58378N64M33114N8M
//     */
//
//    ParserSAM::parse("tests/data/RnaQuin/introns.sam", [&](const Alignment &align, const ParserProgress &)
//    {
//        aligns.push_back(align);
//    });
//
//    REQUIRE(aligns.size() == 5);
//
//    REQUIRE(!aligns[0].spliced);
//    REQUIRE(aligns[0].l.length() == 29);
//    REQUIRE(aligns[0].l == Locus(388598, 388626));
//
//    REQUIRE(aligns[1].spliced);
//    REQUIRE(aligns[1].l.length() == 58378);
//    REQUIRE(aligns[1].l == Locus(388627, 447004));
//
//    REQUIRE(!aligns[2].spliced);
//    REQUIRE(aligns[2].l.length() == 64);
//    REQUIRE(aligns[2].l == Locus(447005, 447068));
//
//    REQUIRE(aligns[3].spliced);
//    REQUIRE(aligns[3].l.length() == 33114);
//    REQUIRE(aligns[3].l == Locus(447069, 480182));
//    
//    REQUIRE(!aligns[4].spliced);
//    REQUIRE(aligns[4].l.length() == 8);
//    REQUIRE(aligns[4].l == Locus(480183, 480190));
//}