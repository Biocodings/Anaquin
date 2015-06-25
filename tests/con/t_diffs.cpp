#include <catch.hpp>
#include "con/c_diffs.hpp"

using namespace Spike;

TEST_CASE("Conjoint_Diffs_Test")
{
    
    CDiffs::analyze("A.sam", "B.sam");

    
    const auto r = CDiffs::analyze("tests/data/con/aligned_A.sam", "tests/data/con/aligned_B.sam");
}