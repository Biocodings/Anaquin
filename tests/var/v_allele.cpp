#include <catch.hpp>
#include "unit/test.hpp"
#include "variant/v_allele.hpp"

using namespace Anaquin;

TEST_CASE("VarAllele_")
{
    Test::variant();

    const auto r  = VAllele::analyze("tests/data/V_1000/VARMXA.approx100xCov.Hg19_with_chrT.given_alleles.SNPs.vcf");
    
    
    
   // const auto lm = r.linear();

   // REQUIRE(r.m.sp()  == Approx(0.977778));
   // REQUIRE(r.m.sn()  == Approx(0.347826));
    //REQUIRE(r.covered == Approx(0.0553359684));
}