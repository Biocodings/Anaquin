#include <catch.hpp>
#include "expression.hpp"

TEST_CASE("D1_Isoforms")
{
    const auto r = Expression::analyze("/Users/tedwong/Sources/QA/Tests/Data/d1/isoforms.fpkm_tracking", ExpressionMode::IsoformExpress);
}

TEST_CASE("D1_Gene")
{
    const auto r = Expression::analyze("/Users/tedwong/Sources/QA/Tests/Data/d1/genes.fpkm_tracking", ExpressionMode::GeneExpress);
}