#include <catch.hpp>
#include "parsers/parser_cdiffs.hpp"

using namespace Spike;

TEST_CASE("Parse_Genes_FPKM_Tracking")
{
    ParserCDiffs::parse("tests/data/rna_sims/gene_exp.diff", [&](const TrackingDiffs &t, const ParserProgress &)
    {
        // Empty Implementation
    });
}

TEST_CASE("Parse_Isoforms_FPKM_Tracking")
{
    ParserCDiffs::parse("tests/data/rna_sims/isoform_exp.diff", [&](const TrackingDiffs &t, const ParserProgress &)
    {
        // Empty Implementation
    });
}