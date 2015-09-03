#include <catch.hpp>
#include "unit/test.hpp"
#include "variant/v_allele.hpp"

using namespace Anaquin;

TEST_CASE("VarAllele_V_1000")
{
    Test::variant();

    const auto r1 = Test::test("-t VarAllele -m data/var/VARMixture_3.0.csv -rbed data/var/VARVariant_1.0.bed -uvcf tests/data/V_1000/VARMXA.approx100xCov.Hg19_with_chrT.given_alleles.SNPs.vcf");

    REQUIRE(r1.status == 0);

    Test::variant();

    const auto r2  = VAllele::analyze("tests/data/V_1000/VARMXA.approx100xCov.Hg19_with_chrT.given_alleles.SNPs.vcf");

   // const auto lm = r.linear();

   // REQUIRE(r.m.sp()  == Approx(0.977778));
   // REQUIRE(r.m.sn()  == Approx(0.347826));
    //REQUIRE(r.covered == Approx(0.0553359684));
}