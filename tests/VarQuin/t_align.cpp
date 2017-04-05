#include <catch.hpp>
#include "test.hpp"
#include "VarQuin/v_align.hpp"

using namespace Anaquin;

TEST_CASE("VAlign_ToyExample")
{
    Test::clear();

    REQUIRE_NOTHROW(Test::test("VarAlign -report 1 -rbed data/VarQuin/A.R.9.bed -ufiles tests/VarQuin_G_NA12878.bam -ufiles tests/VarQuin_G_Sequins.bam"));
}
