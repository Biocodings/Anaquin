#include <catch.hpp>
#include "unit/test.hpp"
#include "trans/t_align.hpp"

using namespace Anaquin;

TEST_CASE("TAlign_All_AllRepeats")
{
    Test::transA();
    
    std::vector<Alignment> aligns;
    
    /*
     * Create synthetic alignments that have mapping only to R2_24
     */
    
    for (auto i = 0; i < 100; i++)
    {
        Alignment align;
        
        align.id      = "chrT";
        align.qName   = "chrT";
        align.i       = 0;
        align.mapped  = true;
        align.spliced = false;
        align.l       = Locus(1122620, 1122629);
        
        aligns.push_back(align);
    }
    
    const auto r  = TAlign::analyze(aligns);
    const auto se = r.data.at(ChrT).eInters.stats();
    const auto si = r.data.at(ChrT).iInters.stats();
    
    REQUIRE(r.data.at(ChrT).unknowns.size() == 0);
    
    REQUIRE(r.data.at(ChrT).overB.h.size()  == 76);
    REQUIRE(r.data.at(ChrT).histE.size() == 76);
    REQUIRE(r.data.at(ChrT).histI.size() == 76);
    
    REQUIRE(se.covered() == Approx(0.0000458388));
    REQUIRE(si.covered() == 0.0);
    
    REQUIRE(r.missing(ChrT, TAlign::Stats::MissingMetrics::MissingExon).i   == 1188);
    REQUIRE(r.missing(ChrT, TAlign::Stats::MissingMetrics::MissingExon).n   == 1190);
    REQUIRE(r.missing(ChrT, TAlign::Stats::MissingMetrics::MissingGene).i   == 76);
    REQUIRE(r.missing(ChrT, TAlign::Stats::MissingMetrics::MissingGene).n   == 76);
    REQUIRE(r.missing(ChrT, TAlign::Stats::MissingMetrics::MissingIntron).i == 1028);
    REQUIRE(r.missing(ChrT, TAlign::Stats::MissingMetrics::MissingIntron).n == 1028);
    
    REQUIRE(r.sn(ChrT, TAlign::Stats::AlignMetrics::AlignBase) == Approx(0.0000504226));
    REQUIRE(r.pc(ChrT, TAlign::Stats::AlignMetrics::AlignBase) == 1.0);
    REQUIRE(r.data.at(ChrT).overB.m.nr() == 218156);
    REQUIRE(r.data.at(ChrT).overB.m.nq() == 10);
    REQUIRE(r.data.at(ChrT).overB.m.tp() == 10);
    REQUIRE(r.data.at(ChrT).overB.m.fp() == 0);
    REQUIRE(r.data.at(ChrT).overB.m.fn() == 218146);
    
    REQUIRE(r.pc(ChrT, TAlign::Stats::AlignMetrics::AlignExon) == 1.0);
    REQUIRE(r.data.at(ChrT).overE.aTP   == 200);
    REQUIRE(r.data.at(ChrT).overE.aFP   == 0);
    REQUIRE(r.data.at(ChrT).overE.aNQ() == 200);
    REQUIRE(r.data.at(ChrT).overE.lTP   == 2);
    REQUIRE(r.data.at(ChrT).overE.lNR   == 1190);
    REQUIRE(r.data.at(ChrT).overE.lFN() == 1188);
    
    REQUIRE(isnan(r.pc(ChrT, TAlign::Stats::AlignMetrics::AlignIntron)));
    REQUIRE(r.data.at(ChrT).overI.aTP   == 0);
    REQUIRE(r.data.at(ChrT).overI.aFP   == 0);
    REQUIRE(r.data.at(ChrT).overI.aNQ() == 0);
    REQUIRE(r.data.at(ChrT).overI.lTP   == 0);
    REQUIRE(r.data.at(ChrT).overI.lNR   == 1028);
    REQUIRE(r.data.at(ChrT).overI.lFN() == 1028);
    
    REQUIRE(r.sn(ChrT, TAlign::Stats::AlignMetrics::AlignExon)   == Approx(0.0016806723));
    REQUIRE(r.sn(ChrT, TAlign::Stats::AlignMetrics::AlignIntron) == 0);
    
    for (auto &i : r.data.at(ChrT).histE)
    {
        if (i.first == "R2_24")
        {
            REQUIRE(i.second == 200);
        }
        else
        {
            REQUIRE(i.second == 0);
        }
    }
    
    for (auto &i : r.data.at(ChrT).histI)
    {
        REQUIRE(i.second == 0);
    }
    
    for (auto &i : r.data.at(ChrT).geneE)
    {
        REQUIRE(i.second.lNR);
        
        if (i.first == "R2_24")
        {
            REQUIRE(r.sn(ChrT, "R2_24")  == Approx(0.0408163265));
            REQUIRE(i.second.pc()  == Approx(1.0));
            REQUIRE(i.second.sn()  == Approx(0.0408163265));
            REQUIRE(i.second.aTP   == 200);
            REQUIRE(i.second.aFP   == 0);
            REQUIRE(i.second.aNQ() == 200);
            REQUIRE(i.second.lTP   == 2);
            REQUIRE(i.second.lNR   == 49);
        }
        else
        {
            REQUIRE(isnan(i.second.pc()));
            REQUIRE(i.second.sn()  == 0);
            REQUIRE(i.second.aTP   == 0);
            REQUIRE(i.second.aFP   == 0);
            REQUIRE(i.second.aNQ() == 0);
            REQUIRE(i.second.lTP   == 0);
        }
    }
    
    for (auto &i : r.data.at(ChrT).geneI)
    {
        if (i.second.lNR)
        {
            REQUIRE(isnan(i.second.pc()));
            REQUIRE(i.second.sn()  == 0);
            REQUIRE(i.second.aTP   == 0);
            REQUIRE(i.second.aFP   == 0);
            REQUIRE(i.second.aNQ() == 0);
            REQUIRE(i.second.lTP   == 0);
        }
        else
        {
            REQUIRE(isnan(i.second.pc()));
            REQUIRE(isnan(i.second.sn()));
            REQUIRE(i.second.aTP   == 0);
            REQUIRE(i.second.aFP   == 0);
            REQUIRE(i.second.aNQ() == 0);
            REQUIRE(i.second.lTP   == 0);
        }
    }
    
    for (auto &i : r.data.at(ChrT).geneB)
    {
        if (i.first == "R2_24")
        {
            REQUIRE(i.second.sn() == Approx(0.0014764506));
            REQUIRE(i.second.pc() == 1.0);
            REQUIRE(i.second.nr() == 6773);
            REQUIRE(i.second.tp() == 10);
            REQUIRE(i.second.fp() == 0);
            REQUIRE(i.second.nq() == 10);
            REQUIRE(i.second.fn() == 6763);
        }
        else
        {
            REQUIRE(i.second.sn() == 0);
            REQUIRE(isnan(i.second.pc()));
            REQUIRE(i.second.nr() != 0);
            REQUIRE(i.second.nq() == 0);
            REQUIRE(i.second.tp() == 0);
            REQUIRE(i.second.fp() == 0);
            REQUIRE(i.second.fn() == i.second.nr());
        }
    }
}

TEST_CASE("TAlign_R2_33_1")
{
    /*
     * R2_33 is a single isoform sequin. We'll generate alignment that covers up the entire sequin.
     * There're two exons in the sequin.
     */
    
    Test::transA();
    
    std::vector<Alignment> aligns;
    
    for (auto i = 0; i < 100; i++)
    {
        Alignment align;
        
        align.id      = "chrT";
        align.qName   = "chrT";
        align.i       = 0;
        align.mapped  = true;
        align.spliced = false;
        
        // The first half covers the first exon while the second half covers the second exon
        align.l = i <= 49 ? Locus(3621204, 3621284) : Locus(3625759, 3625960);
        
        aligns.push_back(align);
    }
    
    const auto r = TAlign::analyze(aligns);
    
    REQUIRE(r.data.at(ChrT).unknowns.size() == 0);
    
    REQUIRE(r.data.at(ChrT).overB.h.size() == 76);
    REQUIRE(r.data.at(ChrT).histE.size()   == 76);
    REQUIRE(r.data.at(ChrT).histI.size()   == 76);
    
    Base sums = 0;
    Base mapped = 0;
    
    for (const auto &i : r.data.at(ChrT).eInters.data())
    {
        if (i.first != "chrT_R2_33_R2_33_1_3621204_3621284" && i.first != "chrT_R2_33_R2_33_1_3625759_3625960")
        {
            REQUIRE(i.second.stats().covered() == 0.00);
        }
        else
        {
            REQUIRE(i.second.stats().covered() == 1.00);
            mapped += i.second.l().length();
        }
        
        sums += i.second.l().length();
    }
    
    const auto covered = static_cast<double>(mapped) / sums;
    
    const auto se = r.data.at(ChrT).eInters.stats();
    const auto si = r.data.at(ChrT).iInters.stats();
    
    REQUIRE(se.covered() == Approx(covered));
    REQUIRE(se.covered() == Approx(0.0012972368));
    REQUIRE(si.covered() == 0.0);
    
    REQUIRE(r.sn(ChrT, TAlign::Stats::AlignMetrics::AlignExon) == Approx(0.0016806723));
    REQUIRE(r.pc(ChrT, TAlign::Stats::AlignMetrics::AlignExon) == 1.0);
    REQUIRE(r.sn(ChrT, TAlign::Stats::AlignMetrics::AlignIntron) == 0);
    REQUIRE(isnan(r.pc(ChrT, TAlign::Stats::AlignMetrics::AlignIntron)));
    REQUIRE(r.sn(ChrT, TAlign::Stats::AlignMetrics::AlignBase) == Approx(0.0012972368));
    REQUIRE(r.pc(ChrT, TAlign::Stats::AlignMetrics::AlignBase) == Approx(1.0));

    REQUIRE(r.data.at(ChrT).overB.m.nr() == 218156);
    REQUIRE(r.data.at(ChrT).overB.m.tp() == mapped);
    REQUIRE(r.data.at(ChrT).overB.m.tp() == 283);
    REQUIRE(r.data.at(ChrT).overB.m.fp() == 0);
    REQUIRE(r.data.at(ChrT).overB.m.fn() == 217873);
    REQUIRE(r.data.at(ChrT).overB.m.nq() == 283);
    
    for (auto &i : r.data.at(ChrT).histI)
    {
        REQUIRE(i.second == 0);
    }
    
    for (auto &i : r.data.at(ChrT).geneE)
    {
        REQUIRE(i.second.lNR);
        
        if (i.first == "R2_33")
        {
            REQUIRE(i.second.pc() == Approx(1.0));
            REQUIRE(i.second.sn()  == Approx(1.0));
            REQUIRE(i.second.aTP   == 100);
            REQUIRE(i.second.aFP   == 0);
            REQUIRE(i.second.aNQ() == 100);
            REQUIRE(i.second.lTP   == 2);
            REQUIRE(i.second.lNR   == 2);
        }
        else
        {
            REQUIRE(isnan(i.second.pc()));
            REQUIRE(i.second.sn()  == 0);
            REQUIRE(i.second.aTP   == 0);
            REQUIRE(i.second.aFP   == 0);
            REQUIRE(i.second.aNQ() == 0);
            REQUIRE(i.second.lTP   == 0);
        }
    }
    
    for (auto &i : r.data.at(ChrT).geneI)
    {
        if (i.second.lNR)
        {
            REQUIRE(isnan(i.second.pc()));
            REQUIRE(i.second.sn()  == 0);
            REQUIRE(i.second.aTP   == 0);
            REQUIRE(i.second.aFP   == 0);
            REQUIRE(i.second.aNQ() == 0);
            REQUIRE(i.second.lTP   == 0);
        }
        else
        {
            REQUIRE(isnan(i.second.pc()));
            REQUIRE(isnan(i.second.sn()));
            REQUIRE(i.second.aTP   == 0);
            REQUIRE(i.second.aFP   == 0);
            REQUIRE(i.second.aNQ() == 0);
            REQUIRE(i.second.lTP   == 0);
        }
    }
    
    for (auto &i : r.data.at(ChrT).geneB)
    {
        if (i.first == "R2_33")
        {
            REQUIRE(i.second.sn() == Approx(1.0));
            REQUIRE(i.second.pc() == 1.0);
            REQUIRE(i.second.nr() == 283);
            REQUIRE(i.second.tp() == 283);
            REQUIRE(i.second.fp() == 0);
            REQUIRE(i.second.nq() == 283);
            REQUIRE(i.second.fn() == 0);
        }
        else
        {
            REQUIRE(i.second.sn() == 0);
            REQUIRE(isnan(i.second.pc()));
            REQUIRE(i.second.nr() != 0);
            REQUIRE(i.second.nq() == 0);
            REQUIRE(i.second.tp() == 0);
            REQUIRE(i.second.fp() == 0);
            REQUIRE(i.second.fn() == i.second.nr());
        }
    }
}

TEST_CASE("TAlign_All_FalsePositives")
{
    Test::transA();
    
    std::vector<Alignment> aligns;
    
    /*
     * Create synthetic alignments that have no mapping to any sequin
     */
    
    for (auto i = 0; i < 100; i++)
    {
        Alignment align;
        
        align.id      = "chrT";
        align.qName   = "chrT";
        align.i       = 0;
        align.mapped  = true;
        align.spliced = false;
        align.l       = Locus(1, 1);
        
        aligns.push_back(align);
    }
    
    const auto r = TAlign::analyze(aligns);
    
    REQUIRE(r.data.at(ChrT).unknowns.size() == 100);
    
    /*
     * There're 76 genes, remember TransAlign does everything at the gene level to avoid
     * the complications due to alternative splicing.
     */
    
    REQUIRE(r.data.at(ChrT).overB.h.size() == 76);
    REQUIRE(r.data.at(ChrT).histE.size() == 76);
    REQUIRE(r.data.at(ChrT).histI.size() == 76);
    
    REQUIRE(r.sn(ChrT, TAlign::Stats::AlignMetrics::AlignExon) == 0);
    REQUIRE(r.pc(ChrT, TAlign::Stats::AlignMetrics::AlignExon) == 0);
    REQUIRE(r.sn(ChrT, TAlign::Stats::AlignMetrics::AlignIntron) == 0);
    REQUIRE(isnan(r.pc(ChrT, TAlign::Stats::AlignMetrics::AlignIntron)));
    
    REQUIRE(r.data.at(ChrT).overB.m.sn() == 0);
    REQUIRE(isnan(r.data.at(ChrT).overB.m.pc()));
    REQUIRE(r.data.at(ChrT).overB.m.nr() == 218156);
    REQUIRE(r.data.at(ChrT).overB.m.nq() == 0);
    REQUIRE(r.data.at(ChrT).overB.m.tp() == 0);
    REQUIRE(r.data.at(ChrT).overB.m.fp() == 0);
    REQUIRE(r.data.at(ChrT).overB.m.fn() == 218156);

    REQUIRE(r.data.at(ChrT).overE.aTP   == 0);
    REQUIRE(r.data.at(ChrT).overE.aFP   == 100);
    REQUIRE(r.data.at(ChrT).overE.aNQ() == 100);
    REQUIRE(r.data.at(ChrT).overE.lTP   == 0);
    REQUIRE(r.data.at(ChrT).overE.lNR   == 1190);
    REQUIRE(r.data.at(ChrT).overE.lFN() == 1190);
    
    REQUIRE(r.data.at(ChrT).overI.aTP   == 0);
    REQUIRE(r.data.at(ChrT).overI.aFP   == 0);
    REQUIRE(r.data.at(ChrT).overI.aNQ() == 0);
    REQUIRE(r.data.at(ChrT).overI.lTP   == 0);
    REQUIRE(r.data.at(ChrT).overI.lNR   == 1028);
    REQUIRE(r.data.at(ChrT).overI.lFN() == 1028);

    REQUIRE(r.sn(ChrT, TAlign::Stats::AlignMetrics::AlignBase) == 0);
    REQUIRE(r.pc(ChrT, TAlign::Stats::AlignMetrics::AlignExon) == 0);
    REQUIRE(r.data.at(ChrT).overB.m.nr() == 218156);
    REQUIRE(r.data.at(ChrT).overB.m.nq() == 0);
    REQUIRE(r.data.at(ChrT).overB.m.tp() == 0);
    REQUIRE(r.data.at(ChrT).overB.m.fp() == 0);
    REQUIRE(r.data.at(ChrT).overB.m.fn() == 218156);
    
    for (auto &i : r.data.at(ChrT).histI)
    {
        REQUIRE(i.second == 0);
    }
    
    for (auto &i : r.data.at(ChrT).geneE)
    {
        REQUIRE(i.second.lNR);
        REQUIRE(isnan(i.second.pc()));
        REQUIRE(i.second.sn()  == 0);
        REQUIRE(i.second.aTP   == 0);
        REQUIRE(i.second.aFP   == 0);
        REQUIRE(i.second.aNQ() == 0);
        REQUIRE(i.second.lTP   == 0);
    }
    
    for (auto &i : r.data.at(ChrT).geneI)
    {
        if (i.second.lNR)
        {
            REQUIRE(isnan(i.second.pc()));
            REQUIRE(i.second.sn()  == 0);
            REQUIRE(i.second.aTP   == 0);
            REQUIRE(i.second.aFP   == 0);
            REQUIRE(i.second.aNQ() == 0);
            REQUIRE(i.second.lTP   == 0);
        }
        else
        {
            REQUIRE(isnan(i.second.pc()));
            REQUIRE(isnan(i.second.sn()));
            REQUIRE(i.second.aTP   == 0);
            REQUIRE(i.second.aFP   == 0);
            REQUIRE(i.second.aNQ() == 0);
            REQUIRE(i.second.lTP   == 0);
        }
    }
    
    for (auto &i : r.data.at(ChrT).geneB)
    {
        REQUIRE(i.second.sn() == 0);
        REQUIRE(isnan(i.second.pc()));
        REQUIRE(i.second.nr() != 0);
        REQUIRE(i.second.nq() == 0);
        REQUIRE(i.second.tp() == 0);
        REQUIRE(i.second.fp() == 0);
        REQUIRE(i.second.fn() == i.second.nr());
    }
}