#include <catch.hpp>
#include "tools/vcf_data.hpp"

using namespace Anaquin;

TEST_CASE("VCF_Synthetic")
{
    const auto r = vcfData(Reader("data/VarQuin/AVA009.v032.vcf"));
    
    REQUIRE(r.countInd()    == 108);
    REQUIRE(r.countIndSyn() == 108);
    REQUIRE(r.countIndGen() == 0);

    REQUIRE(r.countSNP()    == 137);
    REQUIRE(r.countSNPSyn() == 137);
    REQUIRE(r.countSNPGen() == 0);    
}

#ifndef FAST_TEST

/*
 * ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/pilot_data/release/2010_07/trio/snps/CEU.trio.2010_03.genotypes.vcf.gz
 */

TEST_CASE("VCF_CEU_TRIO")
{
    const auto r = vcfData(Reader("data/VarQuin/CEU.trio.2010_03.genotypes.vcf"));
    
    REQUIRE(r.countInd()    == 0);
    REQUIRE(r.countIndSyn() == 0);
    REQUIRE(r.countIndGen() == 0);

    REQUIRE(r.countSNP()    == 3646764);
    REQUIRE(r.countSNPSyn() == 0);
    REQUIRE(r.countSNPGen() == 3646764);
}

#endif