#include <catch.hpp>
#include "test.hpp"

using namespace Anaquin;

TEST_CASE("VarKExpress_1")
{
    REQUIRE_NOTHROW(runTest("VarKExpress -m tests/A.V.11.csv -useqs tests/sample.sf"));
}
