#include <catch.hpp>
#include "stats/analyzer.hpp"

using namespace Spike;

TEST_CASE("RAnalyzer_RNA_Sequin_Counter")
{
    REQUIRE(!RAnalyzer::sequinCounter().empty());
}

TEST_CASE("RAnalyzer_RNA_Exon_Counter")
{
    REQUIRE(!RAnalyzer::exonCounter().empty());
}

TEST_CASE("RAnalyzer_RNA_Intron_Counter")
{
    REQUIRE(!RAnalyzer::intronCounter().empty());
}

TEST_CASE("RAnalyzer_RNA_Gene_Counter")
{
    REQUIRE(!RAnalyzer::geneCounter().empty());
}