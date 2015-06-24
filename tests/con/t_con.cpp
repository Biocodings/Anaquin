#include <catch.hpp>
#include "con/c_join.hpp"

using namespace Spike;

TEST_CASE("CCon_Test")
{
    const auto r = CJoin::analyze("tests/data/con/test.sam");
    
    
    
}

TEST_CASE("CCon_Cufflinks")
{
    const auto r = CJoin::analyze("tests/data/con/cufflinks.sam");

    REQUIRE(r.hist.count("seq1") == 1);
    REQUIRE(r.hist.at("seq1") == 3307);
}