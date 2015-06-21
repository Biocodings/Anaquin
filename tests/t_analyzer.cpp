#include <catch.hpp>
#include "stats/analyzer.hpp"

using namespace Spike;

TEST_CASE("RAnalyzer_RNA_Sequin_Counter")
{
    REQUIRE(!RAnalyzer::sequinCounter().empty());
}