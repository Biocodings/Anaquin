#include <catch.hpp>
#include "assembly.hpp"

TEST_CASE("Generated")
{
    const auto r = Assembly::analyze("/Users/tedwong/Sources/ABCD/transcripts/transcripts.gtf");

    REQUIRE(r.exon.tp == 355);
    REQUIRE(r.exon.fp == 21);
}