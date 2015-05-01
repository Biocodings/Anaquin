#include <catch.hpp>
#include "parsers/parser_sam.hpp"

using namespace Spike;

TEST_CASE("Test_Junction")
{
    std::vector<Alignment> aligns;
    
    /*
     * Cigar: 79M10250N22M
     */
    
    ParserSAM::parse("tests/data/introns.sam", [&](const Alignment &align)
    {
        aligns.push_back(align);
    });

    REQUIRE(aligns.size() == 3);

    REQUIRE(!aligns[0].spliced);
    REQUIRE(aligns[0].l.start == 1);
    REQUIRE(aligns[0].l.end == 79);
    REQUIRE(aligns[0].l.length() == 79);

    REQUIRE(aligns[1].spliced);
    REQUIRE(aligns[1].l.start == 80);
    REQUIRE(aligns[1].l.end == 10329);
    REQUIRE(aligns[1].l.length() == 10250);

    REQUIRE(!aligns[2].spliced);
    REQUIRE(aligns[2].l.start == 10330);
    REQUIRE(aligns[2].l.end == 10351);
    REQUIRE(aligns[2].l.length() == 22);
}