#include <catch.hpp>
#include "tools/gtf_data.hpp"

using namespace Anaquin;

TEST_CASE("GTF_Synthetic")
{
    const auto r = gtfData(Reader("data/RnaQuin/ATR001.v032.gtf"));
    const auto i = r.gIntervals(ChrT);
    
    REQUIRE(i.size() == 78);
    REQUIRE(i.overlap(Locus(6955490, 6955495)));
    REQUIRE(i.contains(Locus(6955490, 6955495)));
    REQUIRE(!i.contains(Locus(6955480, 6955485)));
    REQUIRE(!i.overlap(Locus(6955480, 6955485)));
    
    REQUIRE(r.countGene()     == 78);
    REQUIRE(r.countGeneSyn()  == 78);
    REQUIRE(r.countTrans()    == 164);
    REQUIRE(r.countTransSyn() == 164);
    REQUIRE(r.countUExon()    == 869);
    REQUIRE(r.countUExonSyn() == 869);
    REQUIRE(r.countUIntr()    == 754);
    REQUIRE(r.countUIntrSyn() == 754);

    //REQUIRE(r.il.at(ChrT).at(Locus(6955730, 6960383)) == 1);
    //REQUIRE(r.il.at(ChrT).at(Locus(2227518, 2235700)) == 3);

    REQUIRE(r.ueInters().size() == 1);
    REQUIRE(r.uiInters().size() == 1);

    const auto eIntrs = r.ueInters(ChrT);
    const auto iIntrs = r.uiInters(ChrT);

    REQUIRE(eIntrs.size() == 869); // 1192 for non-unique
    REQUIRE(iIntrs.size() == 754); // 1028 for non-unique

    REQUIRE(eIntrs.contains(Locus(3014253, 3014390)));
    REQUIRE(eIntrs.contains(Locus(3014255, 3014390)));
    REQUIRE(!eIntrs.contains(Locus(3014252, 3014390)));
    REQUIRE(eIntrs.exact(Locus(3014253, 3014390)));
    REQUIRE(!eIntrs.exact(Locus(3014255, 3014390)));
    
    REQUIRE(iIntrs.contains(Locus(3014391, 3014542)));
    REQUIRE(iIntrs.contains(Locus(3014395, 3014540)));
    REQUIRE(!iIntrs.contains(Locus(3014390, 3014540)));
    REQUIRE(iIntrs.exact(Locus(3014391, 3014542)));
    REQUIRE(!iIntrs.exact(Locus(3014395, 3014540)));
}

TEST_CASE("GTF_Merged")
{
    const auto r = gtfData(Reader("data/RnaQuin/merged.gtf"));

    REQUIRE(r.countGeneSyn()  == 78);
    REQUIRE(r.countTransSyn() == 164);
    REQUIRE(r.countUExonSyn() == 869);
    REQUIRE(r.countUIntrSyn() == 754);

    /*
     * cat merged.gtf | grep -v chrT | cut -f3 | grep gene  | wc
     * cat merged.gtf | grep -v chrT | cut -f3 | grep trans | wc
     * cat merged.gtf | grep -v chrT | cut -f3 | grep exon  | wc
     */
    
    REQUIRE(r.countGeneGen()  == 881);
    REQUIRE(r.countTransGen() == 2449);
    REQUIRE(r.countUIntrGen() == 4213);
    REQUIRE(r.countUExonGen() == 6540);
    
    const auto i = r.gIntervals(ChrT);

    REQUIRE(i.size() == 78);
    REQUIRE(i.overlap(Locus(6955490, 6955495)));
    REQUIRE(i.contains(Locus(6955490, 6955495)));
    REQUIRE(!i.contains(Locus(6955480, 6955485)));
    REQUIRE(!i.overlap(Locus(6955480, 6955485)));

    REQUIRE(r.ueInters().size() == 2);
    REQUIRE(r.ueInters(ChrT).size() == 869);     // 1192 for non-unique
    REQUIRE(r.ueInters("chr21").size() == 6540); // 14011 for non-unique

    REQUIRE(r.uiInters().size() == 2);
    REQUIRE(r.uiInters(ChrT).size() == 754);     // 1028 for non-unique
    REQUIRE(r.uiInters("chr21").size() == 4213); // 11553 for non-unique

    //REQUIRE(r.il.at(ChrT).at(Locus(6955730, 6960383)) == 1);
    //REQUIRE(r.il.at(ChrT).at(Locus(2227518, 2235700)) == 3);
}

TEST_CASE("GTF_Sampled")
{
    const auto r = gtfData(Reader("data/tests/RnaQuin/sampled.gtf"));
    REQUIRE(r.countGene() == 598);
}

#ifdef GENCODE_TEST

/*
 *  wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/gencode.v24.annotation.gtf.gz
 *  gunzip gz
 */

TEST_CASE("GTF_GenCode")
{
    const auto r = gtfData(Reader("tests/data/gencode.v24.annotation.gtf"));
    
    REQUIRE(r.countGeneSyn()  == 0);
    REQUIRE(r.countTransSyn() == 0);
    REQUIRE(r.countExonSyn()  == 0);
    REQUIRE(r.countUExonSyn() == 0);
    REQUIRE(r.countIntrSyn()  == 0);
    REQUIRE(r.countUIntrSyn() == 0);
    
    REQUIRE(r.countGene()     == 60554);
    REQUIRE(r.countTrans()    == 199169);
    REQUIRE(r.countGeneGen()  == 60554);
    REQUIRE(r.countTransGen() == 199169);
    REQUIRE(r.countExonGen()  == 1177311);
    REQUIRE(r.countUExonGen() == 570957);
    REQUIRE(r.countIntrGen()  == 976685);
    REQUIRE(r.countUIntrGen() == 347410);
}

#endif