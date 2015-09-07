#include <catch.hpp>
#include "unit/test.hpp"
#include "trans/t_align.hpp"

using namespace Anaquin;

TEST_CASE("TAlign_T_1000_First_100")
{
    Test::trans();
    
    const auto r = Test::test("-t TransAlign -m data/trans/TransMixture_4.1.csv -rgtf data/trans/TransStandard_1.0.gtf -usam abcd.sam");
    
    REQUIRE(r.status == 0);

    

}