#include <catch.hpp>
#include "parser_cdiffs.hpp"

using namespace Spike;

TEST_CASE("Parse_Genes_FPKM_Tracking")
{
    ParserCDiffs::parse("tests/data/rna_sims_2/diffs/gene_exp.diff", [&](const TrackingDiffs &t)
    {
        // Empty Implementation
    });
}

TEST_CASE("Parse_Isoforms_FPKM_Tracking")
{
    ParserCDiffs::parse("tests/data/rna_sims_2/diffs/isoform_exp.diff", [&](const TrackingDiffs &t)
    {
        // Empty Implementation
    });
}