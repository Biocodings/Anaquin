#include <catch.hpp>
#include "unit/test.hpp"
#include "trans/t_align.hpp"

using namespace Anaquin;

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
    
    REQUIRE(r.chrT->unknowns.size() == 0);
    
    REQUIRE(r.chrT->overB.h.size() == 76);
    REQUIRE(r.chrT->histE.size()   == 76);
    REQUIRE(r.chrT->histI.size()   == 76);
    
    Base sums = 0;
    Base mapped = 0;
    
    for (const auto &i : r.chrT->eInters.data())
    {
        if (i.first != "R2_33_R2_33_1_3621204_3621284" && i.first != "R2_33_R2_33_1_3625759_3625960")
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
    
    const auto se = r.chrT->eInters.stats();
    const auto si = r.chrT->iInters.stats();
    
    REQUIRE(se.covered() == Approx(covered));
    REQUIRE(se.covered() == Approx(0.0012972368));
    REQUIRE(si.covered() == 0.0);
    
    REQUIRE(r.sn(TAlign::Stats::AlignMetrics::AlignExon) == Approx(0.0016806723));
    REQUIRE(r.pc(TAlign::Stats::AlignMetrics::AlignExon) == 1.0);
    REQUIRE(r.sn(TAlign::Stats::AlignMetrics::AlignIntron) == 0);
    REQUIRE(isnan(r.pc(TAlign::Stats::AlignMetrics::AlignIntron)));
    REQUIRE(r.sn(TAlign::Stats::AlignMetrics::AlignBase) == Approx(0.0012972368));
    REQUIRE(r.pc(TAlign::Stats::AlignMetrics::AlignBase) == Approx(1.0));

    REQUIRE(r.chrT->overB.m.nr() == 218156);
    REQUIRE(r.chrT->overB.m.tp() == mapped);
    REQUIRE(r.chrT->overB.m.tp() == 283);
    REQUIRE(r.chrT->overB.m.fp() == 0);
    REQUIRE(r.chrT->overB.m.fn() == 217873);
    REQUIRE(r.chrT->overB.m.nq() == 283);
    
    for (auto &i : r.chrT->histI)
    {
        REQUIRE(i.second == 0);
    }
    
    for (auto &i : r.chrT->geneE)
    {
        REQUIRE(i.second.lNR);
        
        if (i.first == "R2_33")
        {
            REQUIRE(i.second.precise() == Approx(1.0));
            REQUIRE(i.second.sn()  == Approx(1.0));
            REQUIRE(i.second.aTP   == 100);
            REQUIRE(i.second.aFP   == 0);
            REQUIRE(i.second.aNQ() == 100);
            REQUIRE(i.second.lTP   == 2);
            REQUIRE(i.second.lNR   == 2);
        }
        else
        {
            REQUIRE(isnan(i.second.precise()));
            REQUIRE(i.second.sn()  == 0);
            REQUIRE(i.second.aTP   == 0);
            REQUIRE(i.second.aFP   == 0);
            REQUIRE(i.second.aNQ() == 0);
            REQUIRE(i.second.lTP   == 0);
        }
    }
    
    for (auto &i : r.chrT->geneI)
    {
        if (i.second.lNR)
        {
            REQUIRE(isnan(i.second.precise()));
            REQUIRE(i.second.sn()  == 0);
            REQUIRE(i.second.aTP   == 0);
            REQUIRE(i.second.aFP   == 0);
            REQUIRE(i.second.aNQ() == 0);
            REQUIRE(i.second.lTP   == 0);
        }
        else
        {
            REQUIRE(isnan(i.second.precise()));
            REQUIRE(isnan(i.second.sn()));
            REQUIRE(i.second.aTP   == 0);
            REQUIRE(i.second.aFP   == 0);
            REQUIRE(i.second.aNQ() == 0);
            REQUIRE(i.second.lTP   == 0);
        }
    }
    
    for (auto &i : r.chrT->geneB)
    {
        if (i.first == "R2_33")
        {
            REQUIRE(i.second.sn() == Approx(1.0));
            REQUIRE(i.second.ac() == 1.0);
            REQUIRE(i.second.nr() == 283);
            REQUIRE(i.second.tp() == 283);
            REQUIRE(i.second.fp() == 0);
            REQUIRE(i.second.nq() == 283);
            REQUIRE(i.second.fn() == 0);
        }
        else
        {
            REQUIRE(i.second.sn() == 0);
            REQUIRE(isnan(i.second.ac()));
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
    
    REQUIRE(r.chrT->unknowns.size() == 100);
    
    /*
     * There're 76 genes, remember TransAlign does everything at the gene level to avoid
     * the complications due to alternative splicing.
     */
    
    REQUIRE(r.chrT->overB.h.size() == 76);
    REQUIRE(r.chrT->histE.size() == 76);
    REQUIRE(r.chrT->histI.size() == 76);
    
    REQUIRE(r.sn(TAlign::Stats::AlignMetrics::AlignExon) == 0);
    REQUIRE(r.pc(TAlign::Stats::AlignMetrics::AlignExon) == 0);
    REQUIRE(r.sn(TAlign::Stats::AlignMetrics::AlignIntron) == 0);
    REQUIRE(isnan(r.pc(TAlign::Stats::AlignMetrics::AlignIntron)));
    
    REQUIRE(r.chrT->overB.m.sn() == 0);
    REQUIRE(isnan(r.chrT->overB.m.ac()));
    REQUIRE(r.chrT->overB.m.nr() == 218156);
    REQUIRE(r.chrT->overB.m.nq() == 0);
    REQUIRE(r.chrT->overB.m.tp() == 0);
    REQUIRE(r.chrT->overB.m.fp() == 0);
    REQUIRE(r.chrT->overB.m.fn() == 218156);

    REQUIRE(r.chrT->overE.aTP   == 0);
    REQUIRE(r.chrT->overE.aFP   == 100);
    REQUIRE(r.chrT->overE.aNQ() == 100);
    REQUIRE(r.chrT->overE.lTP   == 0);
    REQUIRE(r.chrT->overE.lNR   == 1190);
    REQUIRE(r.chrT->overE.lFN() == 1190);
    
    REQUIRE(r.chrT->overI.aTP   == 0);
    REQUIRE(r.chrT->overI.aFP   == 0);
    REQUIRE(r.chrT->overI.aNQ() == 0);
    REQUIRE(r.chrT->overI.lTP   == 0);
    REQUIRE(r.chrT->overI.lNR   == 1028);
    REQUIRE(r.chrT->overI.lFN() == 1028);

    REQUIRE(r.sn(TAlign::Stats::AlignMetrics::AlignBase) == 0);
    REQUIRE(r.pc(TAlign::Stats::AlignMetrics::AlignExon) == 0);
    REQUIRE(r.chrT->overB.m.nr() == 218156);
    REQUIRE(r.chrT->overB.m.nq() == 0);
    REQUIRE(r.chrT->overB.m.tp() == 0);
    REQUIRE(r.chrT->overB.m.fp() == 0);
    REQUIRE(r.chrT->overB.m.fn() == 218156);
    
    for (auto &i : r.chrT->histI)
    {
        REQUIRE(i.second == 0);
    }
    
    for (auto &i : r.chrT->geneE)
    {
        REQUIRE(i.second.lNR);
        REQUIRE(isnan(i.second.precise()));
        REQUIRE(i.second.sn()  == 0);
        REQUIRE(i.second.aTP   == 0);
        REQUIRE(i.second.aFP   == 0);
        REQUIRE(i.second.aNQ() == 0);
        REQUIRE(i.second.lTP   == 0);
    }
    
    for (auto &i : r.chrT->geneI)
    {
        if (i.second.lNR)
        {
            REQUIRE(isnan(i.second.precise()));
            REQUIRE(i.second.sn()  == 0);
            REQUIRE(i.second.aTP   == 0);
            REQUIRE(i.second.aFP   == 0);
            REQUIRE(i.second.aNQ() == 0);
            REQUIRE(i.second.lTP   == 0);
        }
        else
        {
            REQUIRE(isnan(i.second.precise()));
            REQUIRE(isnan(i.second.sn()));
            REQUIRE(i.second.aTP   == 0);
            REQUIRE(i.second.aFP   == 0);
            REQUIRE(i.second.aNQ() == 0);
            REQUIRE(i.second.lTP   == 0);
        }
    }
    
    for (auto &i : r.chrT->geneB)
    {
        REQUIRE(i.second.sn() == 0);
        REQUIRE(isnan(i.second.ac()));
        REQUIRE(i.second.nr() != 0);
        REQUIRE(i.second.nq() == 0);
        REQUIRE(i.second.tp() == 0);
        REQUIRE(i.second.fp() == 0);
        REQUIRE(i.second.fn() == i.second.nr());
    }
}

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
    const auto se = r.chrT->eInters.stats();
    const auto si = r.chrT->iInters.stats();
    
    REQUIRE(r.chrT->unknowns.size() == 0);
    
    REQUIRE(r.chrT->overB.h.size()  == 76);
    REQUIRE(r.chrT->histE.size() == 76);
    REQUIRE(r.chrT->histI.size() == 76);
    
    REQUIRE(se.covered() == Approx(0.0000458388));
    REQUIRE(si.covered() == 0.0);
    
    REQUIRE(r.missing(TAlign::Stats::MissingMetrics::MissingExon).i   == 1188);
    REQUIRE(r.missing(TAlign::Stats::MissingMetrics::MissingExon).n   == 1190);
    REQUIRE(r.missing(TAlign::Stats::MissingMetrics::MissingGene).i   == 76);
    REQUIRE(r.missing(TAlign::Stats::MissingMetrics::MissingGene).n   == 76);
    REQUIRE(r.missing(TAlign::Stats::MissingMetrics::MissingIntron).i == 1028);
    REQUIRE(r.missing(TAlign::Stats::MissingMetrics::MissingIntron).n == 1028);
    
    REQUIRE(r.sn(TAlign::Stats::AlignMetrics::AlignBase) == Approx(0.0000504226));
    REQUIRE(r.pc(TAlign::Stats::AlignMetrics::AlignBase) == 1.0);
    REQUIRE(r.chrT->overB.m.nr() == 218156);
    REQUIRE(r.chrT->overB.m.nq() == 10);
    REQUIRE(r.chrT->overB.m.tp() == 10);
    REQUIRE(r.chrT->overB.m.fp() == 0);
    REQUIRE(r.chrT->overB.m.fn() == 218146);
    
    REQUIRE(r.pc(TAlign::Stats::AlignMetrics::AlignExon) == 1.0);
    REQUIRE(r.chrT->overE.aTP   == 200);
    REQUIRE(r.chrT->overE.aFP   == 0);
    REQUIRE(r.chrT->overE.aNQ() == 200);
    REQUIRE(r.chrT->overE.lTP   == 2);
    REQUIRE(r.chrT->overE.lNR   == 1190);
    REQUIRE(r.chrT->overE.lFN() == 1188);
    
    REQUIRE(isnan(r.pc(TAlign::Stats::AlignMetrics::AlignIntron)));
    REQUIRE(r.chrT->overI.aTP   == 0);
    REQUIRE(r.chrT->overI.aFP   == 0);
    REQUIRE(r.chrT->overI.aNQ() == 0);
    REQUIRE(r.chrT->overI.lTP   == 0);
    REQUIRE(r.chrT->overI.lNR   == 1028);
    REQUIRE(r.chrT->overI.lFN() == 1028);
    
    REQUIRE(r.sn(TAlign::Stats::AlignMetrics::AlignExon)   == Approx(0.0016806723));
    REQUIRE(r.sn(TAlign::Stats::AlignMetrics::AlignIntron) == 0);
    
    for (auto &i : r.chrT->histE)
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
    
    for (auto &i : r.chrT->histI)
    {
        REQUIRE(i.second == 0);
    }
    
    for (auto &i : r.chrT->geneE)
    {
        REQUIRE(i.second.lNR);
        
        if (i.first == "R2_24")
        {
            REQUIRE(r.sn("R2_24")  == Approx(0.0408163265));
            REQUIRE(i.second.precise()  == Approx(1.0));
            REQUIRE(i.second.sn()  == Approx(0.0408163265));
            REQUIRE(i.second.aTP   == 200);
            REQUIRE(i.second.aFP   == 0);
            REQUIRE(i.second.aNQ() == 200);
            REQUIRE(i.second.lTP   == 2);
            REQUIRE(i.second.lNR   == 49);
        }
        else
        {
            REQUIRE(isnan(i.second.precise()));
            REQUIRE(i.second.sn()  == 0);
            REQUIRE(i.second.aTP   == 0);
            REQUIRE(i.second.aFP   == 0);
            REQUIRE(i.second.aNQ() == 0);
            REQUIRE(i.second.lTP   == 0);
        }
    }
    
    for (auto &i : r.chrT->geneI)
    {
        if (i.second.lNR)
        {
            REQUIRE(isnan(i.second.precise()));
            REQUIRE(i.second.sn()  == 0);
            REQUIRE(i.second.aTP   == 0);
            REQUIRE(i.second.aFP   == 0);
            REQUIRE(i.second.aNQ() == 0);
            REQUIRE(i.second.lTP   == 0);
        }
        else
        {
            REQUIRE(isnan(i.second.precise()));
            REQUIRE(isnan(i.second.sn()));
            REQUIRE(i.second.aTP   == 0);
            REQUIRE(i.second.aFP   == 0);
            REQUIRE(i.second.aNQ() == 0);
            REQUIRE(i.second.lTP   == 0);
        }
    }
    
    for (auto &i : r.chrT->geneB)
    {
        if (i.first == "R2_24")
        {
            REQUIRE(i.second.sn() == Approx(0.0014764506));
            REQUIRE(i.second.ac() == 1.0);
            REQUIRE(i.second.nr() == 6773);
            REQUIRE(i.second.tp() == 10);
            REQUIRE(i.second.fp() == 0);
            REQUIRE(i.second.nq() == 10);
            REQUIRE(i.second.fn() == 6763);
        }
        else
        {
            REQUIRE(i.second.sn() == 0);
            REQUIRE(isnan(i.second.ac()));
            REQUIRE(i.second.nr() != 0);
            REQUIRE(i.second.nq() == 0);
            REQUIRE(i.second.tp() == 0);
            REQUIRE(i.second.fp() == 0);
            REQUIRE(i.second.fn() == i.second.nr());
        }
    }
}
