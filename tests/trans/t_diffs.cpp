#include <catch.hpp>
#include "unit/test.hpp"
#include "trans/t_diffs.hpp"

using namespace Anaquin;

TEST_CASE("TDiffs_T_1000_Isoforms")
{
    Test::trans();

//    const auto r1 = Test::test("-t TransDiff -m data/trans/TransMixture_4.1.csv -rgtf data/trans/TransStandard_1.0.gtf -ugdiff tests/data/T_1000/A/G/transcripts.gtf");

    TDiffs::Options o;
    o.level = TDiffs::Isoform;
    
    const auto r  = TDiffs::report("tests/data/T_1000/isoform_exp.diff", o);
    const auto lm = r.linear();
    
    REQUIRE(lm.m  == Approx(0.970367203));
    REQUIRE(lm.r  == Approx(0.8538321793));
    REQUIRE(lm.r2 == Approx(0.7290293904));
}

TEST_CASE("TDiffs_T_1000_Genes")
{
    Test::trans();
    
    TDiffs::Options o;
    o.level = TDiffs::Gene;

    const auto r  = TDiffs::report("tests/data/T_1000/gene_exp.diff", o);
    const auto lm = r.linear();
    
    REQUIRE(lm.m  == Approx(0.9929858312));
    REQUIRE(lm.r  == Approx(0.9951255647));
    REQUIRE(lm.r2 == Approx(0.9902748896));
}