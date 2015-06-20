#include <catch.hpp>
#include "parsers/parser_sam.hpp"

using namespace Spike;

TEST_CASE("Test_Junction")
{
    std::vector<Alignment> aligns;
    
    /*
     * Cigar: 388598 - 29M58378N64M33114N8M
     */

    ParserSAM::parse("tests/data/rna/introns.sam", [&](const Alignment &align, const ParserProgress &)
    {
        aligns.push_back(align);
    });

    REQUIRE(aligns.size() == 5);

    REQUIRE(!aligns[0].spliced);
    REQUIRE(aligns[0].l.length() == 29);
    REQUIRE(aligns[0].l == Locus(388598, 388626));

    REQUIRE(aligns[1].spliced);
    REQUIRE(aligns[1].l.length() == 58378);
    REQUIRE(aligns[1].l == Locus(388627, 447004));

    REQUIRE(!aligns[2].spliced);
    REQUIRE(aligns[2].l.length() == 64);
    REQUIRE(aligns[2].l == Locus(447005, 447068));

    REQUIRE(aligns[3].spliced);
    REQUIRE(aligns[3].l.length() == 33114);
    REQUIRE(aligns[3].l == Locus(447069, 480182));
    
    REQUIRE(!aligns[4].spliced);
    REQUIRE(aligns[4].l.length() == 8);
    REQUIRE(aligns[4].l == Locus(480183, 480190));
}