#include <catch.hpp>
#include "unit/test.hpp"
#include "variant/v_allele.hpp"

using namespace Anaquin;

TEST_CASE("VarAllele_V_1001")
{
    Test::variantA();
    
    const auto r1 = Test::test("-t VarAllele -m data/var/MVA011.v013.csv -rbed data/var/AVA008.v032.bed -uvcf tests/data/V_1001/variant.vcf");

    REQUIRE(r1.status == 0);
    
    Test::variantA();
    
    const auto r2 = VAllele::report("tests/data/V_1001/variant.vcf");
    const auto lm = r2.linear();

    REQUIRE(lm.m  == Approx(0.6937488923));
    REQUIRE(lm.r  == Approx(0.7970460421));
    REQUIRE(lm.r2 == Approx(0.6352823932));
}

TEST_CASE("VarAllele_V_1000")
{
    Test::variantA();

    const auto r1 = Test::test("-t VarAllele -m data/var/MVA011.v013.csv -rbed data/var/AVA008.v032.bed -uvcf tests/data/V_1000/VARMXA.approx100xCov.Hg19_with_chrT.given_alleles.SNPs.vcf");

    REQUIRE(r1.status == 0);

    Test::variantA();

    const auto r2 = VAllele::report("tests/data/V_1000/VARMXA.approx100xCov.Hg19_with_chrT.given_alleles.SNPs.vcf");
    const auto lm = r2.linear();

    REQUIRE(lm.m  == Approx(-0.0475476706));
    REQUIRE(lm.r  == Approx(-0.0630244678));
    REQUIRE(lm.r2 == Approx(0.0039720835));
}