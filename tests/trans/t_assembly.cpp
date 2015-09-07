#include <catch.hpp>
#include "unit/test.hpp"
#include "trans/t_assembly.hpp"

using namespace Anaquin;

TEST_CASE("TAssembly_T_1000")
{
    Test::trans();

    const auto r1 = Test::test("-t TransAssembly -m data/trans/TransMixture_4.1.csv -rgtf data/trans/TransStandard_1.0.gtf -ugtf tests/data/T_1000/A/transcripts.gtf");

    REQUIRE(r1.status == 0);
    
    Test::trans();
    
    TAssembly::Options o;
    
    o.ref   = "data/trans/TransStandard_1.0.gtf";
    o.query = "tests/data/T_1000/A/transcripts.gtf";

    const auto r2 = TAssembly::analyze("tests/data/T_1000/A/transcripts.gtf", o);

    REQUIRE(r2.n_hg38 == 0);
    REQUIRE(r2.n_chrT == 1365);
}