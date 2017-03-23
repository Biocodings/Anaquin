#include <catch.hpp>
#include "test.hpp"
#include "VarQuin/v_align.hpp"

using namespace Anaquin;

TEST_CASE("VAlign_WriteSummary")
{
    VAlign::Stats stats;
    VAlign::Options o;

    REQUIRE_NOTHROW(VAlign::writeSummary("File1", "File2", "File3", stats, o));
    
    o.report = true;
    
    REQUIRE_NOTHROW(VAlign::writeSummary("File1", "File2", "File3", stats, o));
}

TEST_CASE("VAlign_ToyExample")
{
    Test::clear();

    REQUIRE_NOTHROW(Test::test("VarAlign -report 1 -rbed data/VarQuin/A.R.9.bed -ufiles tests/VarQuin_G_NA12878.bam -ufiles tests/VarQuin_G_Sequins.bam"));
}
