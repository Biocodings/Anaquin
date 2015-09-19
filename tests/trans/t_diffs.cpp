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
    
    REQUIRE(lm.m  == Approx(0.3196480888));
    REQUIRE(lm.r  == Approx(0.0313512503));
    REQUIRE(lm.r2 == Approx(0.0009829009));
}

TEST_CASE("TDiffs_T_1000_Genes")
{
    Test::transAB();
    
    TDiffs::Options o;
    o.level = TDiffs::Gene;
    
    const auto r  = TDiffs::report("tests/data/T_1000/gene_exp.diff", o);
    const auto lm = r.linear();
    
    REQUIRE(lm.m  == Approx(0.9691960772));
    REQUIRE(lm.r  == Approx(0.9941174401));
    REQUIRE(lm.r2 == Approx(0.9882694846));
}