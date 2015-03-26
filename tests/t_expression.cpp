#include "gtest/gtest.h"
#include "expression.hpp"

TEST(D1_Isoforms, ExpressionTest)
{
    const auto r = Expression::analyze("/Users/tedwong/Sources/QA/Tests/Data/d1/isoforms.fpkm_tracking", ExpressionMode::IsoformExpress);
}

TEST(D1_Gene, ExpressionTest)
{
    const auto r = Expression::analyze("/Users/tedwong/Sources/QA/Tests/Data/d1/genes.fpkm_tracking", ExpressionMode::GeneExpress);
}