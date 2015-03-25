#include "gtest/gtest.h"
#include "ExpressionAnalyst.hpp"

TEST(D1_Isoforms, ExpressionAnalystTest)
{
    const auto r = ExpressionAnalyst::analyze("/Users/tedwong/Sources/QA/Tests/Data/d1/isoforms.fpkm_tracking", ExpressionMode::IsoformExpress);
}

TEST(D1_Gene, ExpressionAnalystTest)
{
    const auto r = ExpressionAnalyst::analyze("/Users/tedwong/Sources/QA/Tests/Data/d1/genes.fpkm_tracking", ExpressionMode::GeneExpress);
}