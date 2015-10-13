#include <catch.hpp>
#include "variant/v_coverage.hpp"

using namespace Anaquin;

TEST_CASE("VarCoverage_V_1001")
{
    const auto r = VCoverage::stats("tests/data/V_1001/intersect.bam");
}