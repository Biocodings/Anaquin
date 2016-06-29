#include <catch.hpp>
#include "parsers/parser_sam.hpp"

using namespace Anaquin;

TEST_CASE("Test_SoftClip")
{
    std::vector<ParserSAM::Data> r1;
    std::vector<ParserSAM::Info> r2;
    
    ParserSAM::parse("tests/data/test.sam", [&](const ParserSAM::Data &x, const ParserSAM::Info &i)
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