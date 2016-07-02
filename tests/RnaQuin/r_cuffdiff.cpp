#include <catch.hpp>
#include "unit/test.hpp"
#include "RnaQuin/r_cuffdiff.hpp"

using namespace Anaquin;

TEST_CASE("RnaCuffdiff_IsoformExp")
{
    const auto r = RCuffdiff::stats("tests/data/isoform_exp.diff");
    
    REQUIRE(r.data.size() == 17641);

    REQUIRE(r.data[0].cID  == "chr7");
    REQUIRE(r.data[0].gID  == "ENSG00000004059.10");
    REQUIRE(r.data[0].iID  == "ENST00000000233.9");
    REQUIRE(r.data[0].logF == Approx(-0.960368));
    REQUIRE(r.data[0].p    == Approx(0.0019));
    REQUIRE(r.data[0].q    == Approx(0.0152147));
    
    REQUIRE(r.data[1].cID  == "chr12");
    REQUIRE(r.data[1].gID  == "ENSG00000003056.7");
    REQUIRE(r.data[1].iID  == "ENST00000000412.7");
    REQUIRE(r.data[1].logF == Approx(1.30026));
    REQUIRE(r.data[1].p    == Approx(5e-05));
    REQUIRE(r.data[1].q    == Approx(0.00059598));
}

TEST_CASE("RnaCuffdiff_GeneExp")
{
    const auto r = RCuffdiff::stats("tests/data/gene_exp.diff");
    
    REQUIRE(r.data.size()  == 11319);
    
    REQUIRE(r.data[0].cID  == "chr20");
    REQUIRE(r.data[0].gID  == "ENSG00000000419.12");
    REQUIRE(r.data[0].iID  == "");
    REQUIRE(r.data[0].logF == Approx(-0.212994));
    REQUIRE(r.data[0].p    == Approx(0.41945));
    REQUIRE(r.data[0].q    == Approx(0.545719));
    
    REQUIRE(r.data[1].cID  == "chr1");
    REQUIRE(r.data[1].gID  == "ENSG00000000457.13");
    REQUIRE(r.data[1].iID  == "");
    REQUIRE(r.data[1].logF == Approx(1.11264));
    REQUIRE(r.data[1].p    == Approx(0.0179));
    REQUIRE(r.data[1].q    == Approx(0.0463427));
}