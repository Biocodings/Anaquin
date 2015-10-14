#include <catch.hpp>
#include "unit/test.hpp"
#include "trans/t_assembly.hpp"

using namespace Anaquin;

TEST_CASE("TAssembly_T_1000")
{
    Test::transA();

    const auto r1 = Test::test("-t TransAssembly -m data/trans/MTR002.v013.csv -rgtf data/trans/ATR001.v032.gtf -ugtf tests/data/T_1000/A/G/transcripts.gtf");

    REQUIRE(r1.status == 0);
    
    Test::transA();
    
    TAssembly::Options o;
    
    o.ref   = "data/trans/ATR001.v032.gtf";
    o.query = "tests/data/T_1000/A/G/transcripts.gtf";

    const auto r2 = TAssembly::report("tests/data/T_1000/A/G/transcripts.gtf", o);

    REQUIRE(r2.n_expT == 0);
    REQUIRE(r2.n_chrT == 1365);
}