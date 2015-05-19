#include <catch.hpp>
#include "parsers/parser_tmap.hpp"

using namespace Spike;

TEST_CASE("ParserTMap_Cuffcompare")
{
    std::vector<TMap> tmaps;
    
    ParserTMap::parse("tests/data/rna/rna.transcripts.gtf.tmap", [&](const TMap &t, const ParserProgress &)
    {
        tmaps.push_back(t);
    });

    REQUIRE(tmaps.size() == 62);

    REQUIRE(tmaps[0].lFPKM == Approx(2.129797));
    REQUIRE(tmaps[0].fpkm  == Approx(3.122386));
    REQUIRE(tmaps[0].uFPKM == Approx(4.114974));

    REQUIRE(tmaps[1].lFPKM == Approx(10.714247));
    REQUIRE(tmaps[1].fpkm  == Approx(13.773730));
    REQUIRE(tmaps[1].uFPKM == Approx(16.833214));
}