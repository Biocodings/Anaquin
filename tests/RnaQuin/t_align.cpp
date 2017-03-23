//#include <catch.hpp>
//#include "test.hpp"
//#include "RnaQuin/r_align.hpp"
//
//using namespace Anaquin;
//
//typedef RAlign::Stats::AlignMetrics   AlignMetric;
//typedef RAlign::Stats::MissingMetrics MissMetrics;
//
//TEST_CASE("RAlign_All_AllRepeats")
//{
//    Test::transA();
//    std::vector<Alignment> aligns;
//    
//    /*
//     * Create synthetic alignments that have mapping only to R2_24
//     */
//    
//    for (auto i = 0; i < 100; i++)
//    {
//        Alignment align;
//        
//        align.cID     = ChrIS;
//        align.name    = ChrIS;
//        align.i       = 0;
//        align.mapped  = true;
//        align.spliced = false;
//        align.l       = Locus(1122620, 1122629);
//        
//        aligns.push_back(align);
//    }
//    
//    const auto r  = RAlign::analyze(aligns);
//    const auto se = r.data.at(ChrIS).eInters.stats();
//    const auto si = r.data.at(ChrIS).iInters.stats();
//    
//    REQUIRE(r.data.at(ChrIS).unknowns.size() == 0);
//    
//    REQUIRE(r.data.at(ChrIS).overB.hist.size() == 78);
//    REQUIRE(r.data.at(ChrIS).histE.size() == 78);
//    REQUIRE(r.data.at(ChrIS).histI.size() == 78);
//    
//    REQUIRE(se.covered() == Approx(0.0000458388));
//    REQUIRE(si.covered() == 0.0);
//    
//    REQUIRE(r.countMiss(ChrIS, MissMetrics::MissingExon).i   == 1188);
//    REQUIRE(r.countMiss(ChrIS, MissMetrics::MissingExon).n   == 1190);
//    REQUIRE(r.countMiss(ChrIS, MissMetrics::MissingGene).i   == 76);
//    REQUIRE(r.countMiss(ChrIS, MissMetrics::MissingGene).n   == 76);
//    REQUIRE(r.countMiss(ChrIS, MissMetrics::MissingIntron).i == 1028);
//    REQUIRE(r.countMiss(ChrIS, MissMetrics::MissingIntron).n == 1028);
//    
//    REQUIRE(r.sn(ChrIS, AlignMetric::AlignBase) == Approx(0.0000504226));
//    REQUIRE(r.pc(ChrIS, AlignMetric::AlignBase) == 1.0);
//    REQUIRE(r.data.at(ChrIS).overB.m.nr() == 218156);
//    REQUIRE(r.data.at(ChrIS).overB.m.nq() == 10);
//    REQUIRE(r.data.at(ChrIS).overB.m.tp() == 10);
//    REQUIRE(r.data.at(ChrIS).overB.m.fp() == 0);
//    REQUIRE(r.data.at(ChrIS).overB.m.fn() == 218146);
//    
//    REQUIRE(r.pc(ChrIS, RAlign::Stats::AlignMetrics::AlignExon) == 1.0);
//    REQUIRE(r.data.at(ChrIS).overE.aTP   == 200);
//    REQUIRE(r.data.at(ChrIS).overE.aFP   == 0);
//    REQUIRE(r.data.at(ChrIS).overE.aNQ() == 200);
//    REQUIRE(r.data.at(ChrIS).overE.lTP   == 2);
//    REQUIRE(r.data.at(ChrIS).overE.lNR   == 1190);
//    REQUIRE(r.data.at(ChrIS).overE.lFN() == 1188);
//    
//    REQUIRE(isnan(r.pc(ChrIS, AlignMetric::AlignIntron)));
//    REQUIRE(r.data.at(ChrIS).overI.aTP   == 0);
//    REQUIRE(r.data.at(ChrIS).overI.aFP   == 0);
//    REQUIRE(r.data.at(ChrIS).overI.aNQ() == 0);
//    REQUIRE(r.data.at(ChrIS).overI.lTP   == 0);
//    REQUIRE(r.data.at(ChrIS).overI.lNR   == 1028);
//    REQUIRE(r.data.at(ChrIS).overI.lFN() == 1028);
//    
//    REQUIRE(r.sn(ChrIS, AlignMetric::AlignExon)   == Approx(0.0016806723));
//    REQUIRE(r.sn(ChrIS, AlignMetric::AlignIntron) == 0);
//    
//    for (auto &i : r.data.at(ChrIS).histE)
//    {
//        if (i.first == "R2_24")
//        {
//            REQUIRE(i.second == 200);
//        }
//        else
//        {
//            REQUIRE(i.second == 0);
//        }
//    }
//    
//    for (auto &i : r.data.at(ChrIS).histI)
//    {
//        REQUIRE(i.second == 0);
//    }
//    
//    for (auto &i : r.data.at(ChrIS).geneE)
//    {
//        REQUIRE(i.second.lNR);
//        
//        if (i.first == "R2_24")
//        {
//            REQUIRE(r.sn(ChrIS, "R2_24")  == Approx(0.0408163265));
//            REQUIRE(i.second.pc()  == Approx(1.0));
//            REQUIRE(i.second.sn()  == Approx(0.0408163265));
//            REQUIRE(i.second.aTP   == 200);
//            REQUIRE(i.second.aFP   == 0);
//            REQUIRE(i.second.aNQ() == 200);
//            REQUIRE(i.second.lTP   == 2);
//            REQUIRE(i.second.lNR   == 49);
//        }
//        else
//        {
//            REQUIRE(isnan(i.second.pc()));
//            REQUIRE(i.second.sn()  == 0);
//            REQUIRE(i.second.aTP   == 0);
//            REQUIRE(i.second.aFP   == 0);
//            REQUIRE(i.second.aNQ() == 0);
//            REQUIRE(i.second.lTP   == 0);
//        }
//    }
//    
//    for (auto &i : r.data.at(ChrIS).geneI)
//    {
//        if (i.second.lNR)
//        {
//            REQUIRE(isnan(i.second.pc()));
//            REQUIRE(i.second.sn()  == 0);
//            REQUIRE(i.second.aTP   == 0);
//            REQUIRE(i.second.aFP   == 0);
//            REQUIRE(i.second.aNQ() == 0);
//            REQUIRE(i.second.lTP   == 0);
//        }
//        else
//        {
//            REQUIRE(isnan(i.second.pc()));
//            REQUIRE(isnan(i.second.sn()));
//            REQUIRE(i.second.aTP   == 0);
//            REQUIRE(i.second.aFP   == 0);
//            REQUIRE(i.second.aNQ() == 0);
//            REQUIRE(i.second.lTP   == 0);
//        }
//    }
//    
//    for (auto &i : r.data.at(ChrIS).geneB)
//    {
//        if (i.first == "R2_24")
//        {
//            REQUIRE(i.second.sn() == Approx(0.0014764506));
//            REQUIRE(i.second.pc() == 1.0);
//            REQUIRE(i.second.nr() == 6773);
//            REQUIRE(i.second.tp() == 10);
//            REQUIRE(i.second.fp() == 0);
//            REQUIRE(i.second.nq() == 10);
//            REQUIRE(i.second.fn() == 6763);
//        }
//        else
//        {
//            REQUIRE(i.second.sn() == 0);
//            REQUIRE(isnan(i.second.pc()));
//            REQUIRE(i.second.nr() != 0);
//            REQUIRE(i.second.nq() == 0);
//            REQUIRE(i.second.tp() == 0);
//            REQUIRE(i.second.fp() == 0);
//            REQUIRE(i.second.fn() == i.second.nr());
//        }
//    }
//}
//
//TEST_CASE("RAlign_R2_33_1")
//{
//    /*
//     * R2_33 is a single isoform sequin. We'll generate alignment that covers up the entire sequin.
//     * There're two exons in the sequin.
//     */
//    
//    Test::transA();
//    
//    std::vector<Alignment> aligns;
//    
//    for (auto i = 0; i < 100; i++)
//    {
//        Alignment align;
//        
//        align.cID     = ChrIS;
//        align.name    = ChrIS;
//        align.i       = 0;
//        align.mapped  = true;
//        align.spliced = false;
//        
//        // The first half covers the first exon while the second half covers the second exon
//        align.l = i <= 49 ? Locus(3621204, 3621284) : Locus(3625759, 3625960);
//        
//        aligns.push_back(align);
//    }
//    
//    const auto r = RAlign::analyze(aligns);
//    
//    REQUIRE(r.data.at(ChrIS).unknowns.size() == 0);
//    
//    REQUIRE(r.data.at(ChrIS).overB.hist.size() == 76);
//    REQUIRE(r.data.at(ChrIS).histE.size()   == 76);
//    REQUIRE(r.data.at(ChrIS).histI.size()   == 76);
//    
//    Base sums = 0;
//    Base mapped = 0;
//    
//    for (const auto &i : r.data.at(ChrIS).eInters.data())
//    {
//        if (i.first != "ChrIS_R2_33_R2_33_1_3621204_3621284" && i.first != "ChrIS_R2_33_R2_33_1_3625759_3625960")
//        {
//            REQUIRE(i.second.stats().covered() == 0.00);
//        }
//        else
//        {
//            REQUIRE(i.second.stats().covered() == 1.00);
//            mapped += i.second.l().length();
//        }
//        
//        sums += i.second.l().length();
//    }
//    
//    const auto covered = static_cast<double>(mapped) / sums;
//    
//    const auto se = r.data.at(ChrIS).eInters.stats();
//    const auto si = r.data.at(ChrIS).iInters.stats();
//    
//    REQUIRE(se.covered() == Approx(covered));
//    REQUIRE(se.covered() == Approx(0.0012972368));
//    REQUIRE(si.covered() == 0.0);
//    
//    REQUIRE(r.sn(ChrIS, AlignMetric::AlignExon) == Approx(0.0016806723));
//    REQUIRE(r.pc(ChrIS, AlignMetric::AlignExon) == 1.0);
//    REQUIRE(r.sn(ChrIS, AlignMetric::AlignIntron) == 0);
//    REQUIRE(isnan(r.pc(ChrIS, AlignMetric::AlignIntron)));
//    REQUIRE(r.sn(ChrIS, AlignMetric::AlignBase) == Approx(0.0012972368));
//    REQUIRE(r.pc(ChrIS, AlignMetric::AlignBase) == Approx(1.0));
//
//    REQUIRE(r.data.at(ChrIS).overB.m.nr() == 218156);
//    REQUIRE(r.data.at(ChrIS).overB.m.tp() == mapped);
//    REQUIRE(r.data.at(ChrIS).overB.m.tp() == 283);
//    REQUIRE(r.data.at(ChrIS).overB.m.fp() == 0);
//    REQUIRE(r.data.at(ChrIS).overB.m.fn() == 217873);
//    REQUIRE(r.data.at(ChrIS).overB.m.nq() == 283);
//    
//    for (auto &i : r.data.at(ChrIS).histI)
//    {
//        REQUIRE(i.second == 0);
//    }
//    
//    for (auto &i : r.data.at(ChrIS).geneE)
//    {
//        REQUIRE(i.second.lNR);
//        
//        if (i.first == "R2_33")
//        {
//            REQUIRE(i.second.pc() == Approx(1.0));
//            REQUIRE(i.second.sn()  == Approx(1.0));
//            REQUIRE(i.second.aTP   == 100);
//            REQUIRE(i.second.aFP   == 0);
//            REQUIRE(i.second.aNQ() == 100);
//            REQUIRE(i.second.lTP   == 2);
//            REQUIRE(i.second.lNR   == 2);
//        }
//        else
//        {
//            REQUIRE(isnan(i.second.pc()));
//            REQUIRE(i.second.sn()  == 0);
//            REQUIRE(i.second.aTP   == 0);
//            REQUIRE(i.second.aFP   == 0);
//            REQUIRE(i.second.aNQ() == 0);
//            REQUIRE(i.second.lTP   == 0);
//        }
//    }
//    
//    for (auto &i : r.data.at(ChrIS).geneI)
//    {
//        if (i.second.lNR)
//        {
//            REQUIRE(isnan(i.second.pc()));
//            REQUIRE(i.second.sn()  == 0);
//            REQUIRE(i.second.aTP   == 0);
//            REQUIRE(i.second.aFP   == 0);
//            REQUIRE(i.second.aNQ() == 0);
//            REQUIRE(i.second.lTP   == 0);
//        }
//        else
//        {
//            REQUIRE(isnan(i.second.pc()));
//            REQUIRE(isnan(i.second.sn()));
//            REQUIRE(i.second.aTP   == 0);
//            REQUIRE(i.second.aFP   == 0);
//            REQUIRE(i.second.aNQ() == 0);
//            REQUIRE(i.second.lTP   == 0);
//        }
//    }
//    
//    for (auto &i : r.data.at(ChrIS).geneB)
//    {
//        if (i.first == "R2_33")
//        {
//            REQUIRE(i.second.sn() == Approx(1.0));
//            REQUIRE(i.second.pc() == 1.0);
//            REQUIRE(i.second.nr() == 283);
//            REQUIRE(i.second.tp() == 283);
//            REQUIRE(i.second.fp() == 0);
//            REQUIRE(i.second.nq() == 283);
//            REQUIRE(i.second.fn() == 0);
//        }
//        else
//        {
//            REQUIRE(i.second.sn() == 0);
//            REQUIRE(isnan(i.second.pc()));
//            REQUIRE(i.second.nr() != 0);
//            REQUIRE(i.second.nq() == 0);
//            REQUIRE(i.second.tp() == 0);
//            REQUIRE(i.second.fp() == 0);
//            REQUIRE(i.second.fn() == i.second.nr());
//        }
//    }
//}
//
//TEST_CASE("RAlign_All_FalsePositives")
//{
//    Test::transA();
//    
//    std::vector<Alignment> aligns;
//    
//    /*
//     * Create synthetic alignments that have no mapping to any sequin
//     */
//    
//    for (auto i = 0; i < 100; i++)
//    {
//        Alignment align;
//        
//        align.cID     = ChrIS;
//        align.name    = ChrIS;
//        align.i       = 0;
//        align.mapped  = true;
//        align.spliced = false;
//        align.l       = Locus(1, 1);
//
//        aligns.push_back(align);
//    }
//    
//    const auto r = RAlign::analyze(aligns);
//    
//    REQUIRE(r.data.at(ChrIS).unknowns.size() == 100);
//    
//    /*
//     * There're 76 genes, remember RnaAlign does everything at the gene level to avoid
//     * the complications due to alternative splicing.
//     */
//    
//    REQUIRE(r.data.at(ChrIS).overB.hist.size() == 76);
//    REQUIRE(r.data.at(ChrIS).histE.size() == 76);
//    REQUIRE(r.data.at(ChrIS).histI.size() == 76);
//    
//    REQUIRE(r.sn(ChrIS, AlignMetric::AlignExon) == 0);
//    REQUIRE(r.pc(ChrIS, AlignMetric::AlignExon) == 0);
//    REQUIRE(r.sn(ChrIS, AlignMetric::AlignIntron) == 0);
//    REQUIRE(isnan(r.pc(ChrIS, AlignMetric::AlignIntron)));
//    
//    REQUIRE(r.data.at(ChrIS).overB.m.sn() == 0);
//    REQUIRE(isnan(r.data.at(ChrIS).overB.m.pc()));
//    REQUIRE(r.data.at(ChrIS).overB.m.nr() == 218156);
//    REQUIRE(r.data.at(ChrIS).overB.m.nq() == 0);
//    REQUIRE(r.data.at(ChrIS).overB.m.tp() == 0);
//    REQUIRE(r.data.at(ChrIS).overB.m.fp() == 0);
//    REQUIRE(r.data.at(ChrIS).overB.m.fn() == 218156);
//
//    REQUIRE(r.data.at(ChrIS).overE.aTP   == 0);
//    REQUIRE(r.data.at(ChrIS).overE.aFP   == 100);
//    REQUIRE(r.data.at(ChrIS).overE.aNQ() == 100);
//    REQUIRE(r.data.at(ChrIS).overE.lTP   == 0);
//    REQUIRE(r.data.at(ChrIS).overE.lNR   == 1190);
//    REQUIRE(r.data.at(ChrIS).overE.lFN() == 1190);
//    
//    REQUIRE(r.data.at(ChrIS).overI.aTP   == 0);
//    REQUIRE(r.data.at(ChrIS).overI.aFP   == 0);
//    REQUIRE(r.data.at(ChrIS).overI.aNQ() == 0);
//    REQUIRE(r.data.at(ChrIS).overI.lTP   == 0);
//    REQUIRE(r.data.at(ChrIS).overI.lNR   == 1028);
//    REQUIRE(r.data.at(ChrIS).overI.lFN() == 1028);
//
//    REQUIRE(r.sn(ChrIS, AlignMetric::AlignBase) == 0);
//    REQUIRE(r.pc(ChrIS, AlignMetric::AlignExon) == 0);
//    REQUIRE(r.data.at(ChrIS).overB.m.nr() == 218156);
//    REQUIRE(r.data.at(ChrIS).overB.m.nq() == 0);
//    REQUIRE(r.data.at(ChrIS).overB.m.tp() == 0);
//    REQUIRE(r.data.at(ChrIS).overB.m.fp() == 0);
//    REQUIRE(r.data.at(ChrIS).overB.m.fn() == 218156);
//    
//    for (auto &i : r.data.at(ChrIS).histI)
//    {
//        REQUIRE(i.second == 0);
//    }
//    
//    for (auto &i : r.data.at(ChrIS).geneE)
//    {
//        REQUIRE(i.second.lNR);
//        REQUIRE(isnan(i.second.pc()));
//        REQUIRE(i.second.sn()  == 0);
//        REQUIRE(i.second.aTP   == 0);
//        REQUIRE(i.second.aFP   == 0);
//        REQUIRE(i.second.aNQ() == 0);
//        REQUIRE(i.second.lTP   == 0);
//    }
//    
//    for (auto &i : r.data.at(ChrIS).geneI)
//    {
//        if (i.second.lNR)
//        {
//            REQUIRE(isnan(i.second.pc()));
//            REQUIRE(i.second.sn()  == 0);
//            REQUIRE(i.second.aTP   == 0);
//            REQUIRE(i.second.aFP   == 0);
//            REQUIRE(i.second.aNQ() == 0);
//            REQUIRE(i.second.lTP   == 0);
//        }
//        else
//        {
//            REQUIRE(isnan(i.second.pc()));
//            REQUIRE(isnan(i.second.sn()));
//            REQUIRE(i.second.aTP   == 0);
//            REQUIRE(i.second.aFP   == 0);
//            REQUIRE(i.second.aNQ() == 0);
//            REQUIRE(i.second.lTP   == 0);
//        }
//    }
//    
//    for (auto &i : r.data.at(ChrIS).geneB)
//    {
//        REQUIRE(i.second.sn() == 0);
//        REQUIRE(isnan(i.second.pc()));
//        REQUIRE(i.second.nr() != 0);
//        REQUIRE(i.second.nq() == 0);
//        REQUIRE(i.second.tp() == 0);
//        REQUIRE(i.second.fp() == 0);
//        REQUIRE(i.second.fn() == i.second.nr());
//    }
//}