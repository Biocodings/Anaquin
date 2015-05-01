//#include <catch.hpp>
//#include "rna/r_align.hpp"
//
//using namespace Spike;
//
//static std::string exts[] = { "sam", "bam" };
//
//TEST_CASE("RAlign_Simulations_Base")
//{
//    for (auto ex : exts)
//    {
//        const auto r = RAlign::analyze("tests/data/rna_sims/accepted_hits." + ex);
//
//        //REQUIRE(r.n  == 99066);
//        //REQUIRE(r.nr == 99066);
//        //REQUIRE(r.nq == 0);
//
//        REQUIRE(r.sb.id == "R_5_2");
//        REQUIRE(r.sb.counts == 12);
//        REQUIRE(r.sb.abund == 9765.0);
//        
//        REQUIRE(r.se.id == "R_5_2");
//        REQUIRE(r.se.counts == 2);
//        REQUIRE(r.se.abund == 9765.0);
//    }
//}
//
//TEST_CASE("RAlign_Simulations_Exon")
//{
//    for (auto ex : exts)
//    {
//        RAlign::Options options;
//        //options.level = RAlign::Exon;
//        const auto r = RAlign::analyze("tests/data/rna_sims/accepted_hits." + ex, options);
//    
//       // REQUIRE(r.n  == 53394);
//        //REQUIRE(r.nr == 53394);
//       // REQUIRE(r.nq == 0);
//        REQUIRE(r.mb.tp() == 52835);
//        REQUIRE(r.mb.fp() == 559);
//    }
//}
//
//TEST_CASE("RAlign_RNA_Cufflinks")
//{
//    // The sample file was taken from Cufflink's source distribution. It's obviously independent.
//    const auto r = RAlign::analyze("tests/data/cufflinks.sam");
//    
//    //REQUIRE(0 == r.nr);
//   // REQUIRE(3271 == r.n);
//    //REQUIRE(3271 == r.nq);
//    //REQUIRE(0 == r.mb.tp());
//    //REQUIRE(0 == r.mb.fp());
//    //REQUIRE(isnan(r.mb.sp()));
//    //REQUIRE(isnan(r.mb.sn()));
//    //REQUIRE(0 == r.dilution());
//}