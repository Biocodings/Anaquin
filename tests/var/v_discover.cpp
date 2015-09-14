#include <catch.hpp>
#include "unit/test.hpp"
#include "variant/v_discover.hpp"

using namespace Anaquin;

TEST_CASE("VarDiscover_V_1000")
{
    Test::variant();
    
    const auto r = VDiscover::report("tests/data/V_1000/VARMXA.approx100xCov.Hg19_with_chrT.given_alleles.SNPs.vcf");

    REQUIRE(r.m.sp() == Approx(1.0));
    REQUIRE(r.m.sn() == Approx(0.5415019763));
}