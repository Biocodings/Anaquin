#include <catch.hpp>
#include "conj/c_join.hpp"

using namespace Spike;

TEST_CASE("CJoin_RNA_Sequin_Counter")
{
    const auto r = CJoin::analyze("tests/data/conj/cufflinks.sam");

    REQUIRE(r.hist.count("seq1") == 1);
    REQUIRE(r.hist.at("seq1") == 3307);
}