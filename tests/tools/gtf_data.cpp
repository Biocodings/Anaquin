#include <catch.hpp>
#include "tools/gtf_data.hpp"

using namespace Anaquin;

TEST_CASE("GTF_Sampled")
{
    const auto r = gtfData(Reader("data/tests/RnaQuin/sampled.gtf"));

    REQUIRE(r.countGene() == 598);
}

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
    REQUIRE(r.countExon()     == 1192);
    REQUIRE(r.countExonSyn()  == 1192);
    REQUIRE(r.countUExon()    == 869);
    REQUIRE(r.countUExonSyn() == 869);
    REQUIRE(r.countIntr()     == 1028);
    REQUIRE(r.countIntrSyn()  == 1028);
    REQUIRE(r.countUIntr()    == 754);
    REQUIRE(r.countUIntrSyn() == 754);

    REQUIRE(r.il.at(ChrT).at(Locus(6955730, 6960383)) == 1);
    REQUIRE(r.il.at(ChrT).at(Locus(2227518, 2235700)) == 3);

    REQUIRE(r.eIntervals().size() == 1);
    REQUIRE(r.iIntervals().size() == 1);

    const auto eIntrs = r.eIntervals(ChrT);
    const auto iIntrs = r.iIntervals(ChrT);

    REQUIRE(eIntrs.size() == 1192);
    REQUIRE(iIntrs.size() == 1028);

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
    REQUIRE(r.countExonSyn()  == 1192);
    REQUIRE(r.countUExonSyn() == 869);
    REQUIRE(r.countIntrSyn()  == 1028);
    REQUIRE(r.countUIntrSyn() == 754);

    /*
     * cat merged.gtf | grep -v chrT | cut -f3 | grep gene  | wc
     * cat merged.gtf | grep -v chrT | cut -f3 | grep trans | wc
     * cat merged.gtf | grep -v chrT | cut -f3 | grep exon  | wc
     */
    
    REQUIRE(r.countGeneGen()  == 881);
    REQUIRE(r.countTransGen() == 2449);
    REQUIRE(r.countExonGen()  == 14011);
    REQUIRE(r.countIntrGen()  == 11553);
    REQUIRE(r.countUIntrGen() == 4213);
    REQUIRE(r.countUExonGen() == 6540);
    
    const auto i = r.gIntervals(ChrT);

    REQUIRE(i.size() == 78);
    REQUIRE(i.overlap(Locus(6955490, 6955495)));
    REQUIRE(i.contains(Locus(6955490, 6955495)));
    REQUIRE(!i.contains(Locus(6955480, 6955485)));
    REQUIRE(!i.overlap(Locus(6955480, 6955485)));

    REQUIRE(r.eIntervals().size() == 2);
    REQUIRE(r.eIntervals(ChrT).size() == 1192);
    REQUIRE(r.eIntervals("chr21").size() == 14011);

    REQUIRE(r.iIntervals().size() == 2);
    REQUIRE(r.iIntervals(ChrT).size() == 1028);
    REQUIRE(r.iIntervals("chr21").size() == 11553);

    REQUIRE(r.il.at(ChrT).at(Locus(6955730, 6960383)) == 1);
    REQUIRE(r.il.at(ChrT).at(Locus(2227518, 2235700)) == 3);
}

//TEST_CASE("GTF_Summary_2")
//{
//    const auto r = gtfData(Reader("data/tests/RnaQuin/combined.gtf"));
//    
//    REQUIRE(r.countGene()     == 958);
//    REQUIRE(r.countGeneSyn()  == 78);
//    REQUIRE(r.countTrans()    == 2599);
//    REQUIRE(r.countTransSyn() == 164);
//    REQUIRE(r.countExon()     == 15110);
//    REQUIRE(r.countExonSyn()  == 1192);
//    REQUIRE(r.countIntr()     == 12511);
//    REQUIRE(r.countIntrSyn()  == 1028);
//}

/*
 *  wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/gencode.v24.annotation.gtf.gz
 *  gunzip gz
 */

//TEST_CASE("GTF_GenCode")
//{
//    const auto r = gtfData(Reader("tests/data/gencode.v24.annotation.gtf"));
//    
//    REQUIRE(r.countGeneSyn()  == 0);
//    REQUIRE(r.countTransSyn() == 0);
//    REQUIRE(r.countExonSyn()  == 0);
//    REQUIRE(r.countUExonSyn() == 0);
//    REQUIRE(r.countIntrSyn()  == 0);
//    REQUIRE(r.countUIntrSyn() == 0);
//    
//    REQUIRE(r.countGene()     == 60554);
//    REQUIRE(r.countTrans()    == 199169);
//    REQUIRE(r.countGeneGen()  == 60554);
//    REQUIRE(r.countTransGen() == 199169);
//    REQUIRE(r.countExonGen()  == 1177311);
//    REQUIRE(r.countUExonGen() == 570957);
//    REQUIRE(r.countIntrGen()  == 976685);
//    REQUIRE(r.countUIntrGen() == 347410);
//}