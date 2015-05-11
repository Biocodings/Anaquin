#include <catch.hpp>
#include "dna/d_align.hpp"

using namespace Spike;

static std::string exts[] = { "sam", "bam" };

TEST_CASE("DAlign_Simulations")
{
    const auto r = DAlign::analyze("tests/data/dna_sims/accepted_hits.sam");
}