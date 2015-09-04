#include <catch.hpp>
#include "unit/test.hpp"
#include "trans/t_assembly.hpp"

using namespace Anaquin;

#ifdef REGRESSION

TEST_CASE("TAssembly_T_1000")
{
    Test::trans();

    const auto r = Test::test("-t TransAssembly -m data/trans/TransMixture_4.1.csv -rgtf data/trans/TransStandard_1.0.gtf -ugtf /share/Projects/Tutorials/TransQuin/B/G/transcripts.gtf");

    REQUIRE(r.status == 0);
}

#endif
