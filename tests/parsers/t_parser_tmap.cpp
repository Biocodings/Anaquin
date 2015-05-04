#include <catch.hpp>
#include "parsers/parser_tmap.hpp"

using namespace Spike;

TEST_CASE("ParserTMap_Cuffcompare")
{
    std::vector<TMap> tmaps;
    
    ParserTMap::parse("tests/data/rna_sims/rna.transcripts_an.gtf.tmap", [&](const TMap &t, const ParserProgress &)
    {
        tmaps.push_back(t);
    });

    REQUIRE(tmaps.size() == 61);

    REQUIRE(tmaps[0].lFPKM == Approx(48.175752));
    REQUIRE(tmaps[0].fpkm  == Approx(115.882484));
    REQUIRE(tmaps[0].uFPKM == Approx(183.589216));

    REQUIRE(tmaps[1].lFPKM == 0.00);
    REQUIRE(tmaps[1].fpkm  == 0.00);
    REQUIRE(tmaps[1].uFPKM == 0.00);

    REQUIRE(tmaps[5].lFPKM == Approx(2682.960189));
    REQUIRE(tmaps[5].fpkm  == Approx(3414.686897));
    REQUIRE(tmaps[5].uFPKM == Approx(4146.413606));
}