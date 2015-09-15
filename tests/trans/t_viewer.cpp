#include "unit/test.hpp"
#include "trans/t_viewer.hpp"

using namespace Anaquin;

//#ifdef REGRESSION

TEST_CASE("TViewer_T_1001")
{
    const auto r1 = Test::test("-t TransIGV -o ABCD -ubam accepted_hits.bam");
    
    REQUIRE(r1.status == 0);
}

//#endif