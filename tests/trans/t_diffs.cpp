#include <catch.hpp>
#include "unit/test.hpp"
#include "trans/t_diffs.hpp"

using namespace Anaquin;

TEST_CASE("TDiffs_T_1000_Isoforms")
{
    Test::transAB();
    
    TDiffs::Options o;
    o.level = TDiffs::Isoform;
    
    const auto r  = TDiffs::report("tests/data/T_1000/isoform_exp.diff", o);
    const auto lm = r.linear();
    
    REQUIRE(lm.m  == Approx(0.9702549412));
    REQUIRE(lm.r  == Approx(0.8546476527));
    REQUIRE(lm.r2 == Approx(0.7304226103));
}

TEST_CASE("TDiffs_T_1000_Genes")
{
    Test::transAB();
    
    TDiffs::Options o;
    o.level = TDiffs::Gene;
    
    const auto r  = TDiffs::report("tests/data/T_1000/gene_exp.diff", o);
    const auto lm = r.linear();
    
    REQUIRE(lm.m  == Approx(0.9921460476));
    REQUIRE(lm.r  == Approx(0.9953198807));
    REQUIRE(lm.r2 == Approx(0.9906616648));
}