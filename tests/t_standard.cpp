#include <catch.hpp>
#include "data/standard.hpp"

using namespace Spike;

TEST_CASE("Standard_RNA_Sequins")
{
    const auto &s = Standard::instance();
    
    REQUIRE(s.r_seqs_A.count("R1_92_1")      == 1);
    REQUIRE(s.r_seqs_A.at("R1_92_1").l.start == 144926);
    REQUIRE(s.r_seqs_A.at("R1_92_1").l.end   == 2199400);
}

TEST_CASE("Standard_ID")
{
	REQUIRE("chrT" == Standard::instance().id);
}
