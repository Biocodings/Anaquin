#include <catch.hpp>
#include "unit/test.hpp"
#include "trans/t_align.hpp"

using namespace Anaquin;

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
    
    const auto r = TAlign::stats(aligns);
    
    REQUIRE(r.unknowns.size() == 100);
    
    /*
     * There're 76 genes, remember TransAlign does everything at the gene level to avoid
     * the complications due to alternative splicing.
     */
    
    REQUIRE(r.pb.h.size() == 76);
    REQUIRE(r.pe.h.size() == 76);
    REQUIRE(r.pi.h.size() == 76);
    
    REQUIRE(r.pe.m.sn() == 0);
    REQUIRE(r.pe.m.ac() == 0);
    REQUIRE(r.pe.m.nr() == 76);
    REQUIRE(r.pe.m.tp() == 0);
    REQUIRE(r.pe.m.fp() == 100);
    REQUIRE(r.pe.m.fn() == 76);
    
    REQUIRE(r.pi.m.sn() == 0);
    REQUIRE(isnan(r.pi.m.ac()));
    REQUIRE(r.pi.m.nr() == 76);
    REQUIRE(r.pi.m.nq() == 0);
    REQUIRE(r.pi.m.tp() == 0);
    REQUIRE(r.pi.m.fp() == 0);
    REQUIRE(r.pi.m.fn() == 76);
    
    REQUIRE(r.pb.m.sn() == 0);
    REQUIRE(r.pb.m.ac() == 0);
    REQUIRE(r.pb.m.nr() == 218156);
    REQUIRE(r.pb.m.nq() == 1);
    REQUIRE(r.pb.m.tp() == 0);
    REQUIRE(r.pb.m.fp() == 1);
    REQUIRE(r.pb.m.fn() == 218156);
    
    for (auto &i : r.se)
    {
        REQUIRE(i.second.sn() == 0);
        REQUIRE(isnan(i.second.ac()));
        REQUIRE(i.second.nr() != 0); // Number of exons
        REQUIRE(i.second.nq() == 0);
        REQUIRE(i.second.tp() == 0);
        REQUIRE(i.second.fp() == 0);
        REQUIRE(i.second.fn() == i.second.nr());
    }
    
    for (auto &i : r.si)
    {
        // Not all sequin has an intron..
        if (i.second.nr())
        {
            REQUIRE(i.second.sn() == 0);
            REQUIRE(isnan(i.second.ac()));
            REQUIRE(i.second.nr() != 0); // Number of introns
            REQUIRE(i.second.nq() == 0);
            REQUIRE(i.second.tp() == 0);
            REQUIRE(i.second.fp() == 0);
            REQUIRE(i.second.fn() == i.second.nr());
        }
    }
    
    for (auto &i : r.sb)
    {
        REQUIRE(i.second.sn() == 0);
        REQUIRE(isnan(i.second.ac()));
        REQUIRE(i.second.nr() != 0); // Number of bases
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
    
    const auto r  = TAlign::stats(aligns);
    const auto se = r.eInters.stats();
    const auto si = r.iInters.stats();

    REQUIRE(r.unknowns.size() == 0);
  
    REQUIRE(r.pb.h.size() == 76);
    REQUIRE(r.pe.h.size() == 76);
    REQUIRE(r.pi.h.size() == 76);
  
    REQUIRE(se.covered() == Approx(0.0000458388));
    REQUIRE(si.covered() == 0.0);

    REQUIRE(r.pb.m.sn() == Approx(0.0000504226));
    REQUIRE(r.pb.m.ac() == 1.0);
    REQUIRE(r.pb.m.nr() == 218156);
    REQUIRE(r.pb.m.nq() == 10);
    REQUIRE(r.pb.m.tp() == 10);
    REQUIRE(r.pb.m.fp() == 0);
    REQUIRE(r.pb.m.fn() == 218146);

    REQUIRE(r.pe.m.sn() == Approx(0.0131578947));
    REQUIRE(r.pe.m.ac() == 1.0);
    REQUIRE(r.pe.m.nr() == 76);
    REQUIRE(r.pe.m.tp() == 1);
    REQUIRE(r.pe.m.fn() == 75);
    
    REQUIRE(r.pi.m.sn() == 0);
    REQUIRE(isnan(r.pi.m.ac()));
    REQUIRE(r.pi.m.nr() == 76);
    REQUIRE(r.pi.m.tp() == 0);
    REQUIRE(r.pi.m.fn() == 76);
    
    for (auto &i : r.si)
    {
        // Not all sequin has an intron..
        if (i.second.nr())
        {
            REQUIRE(i.second.sn() == 0);
            REQUIRE(isnan(i.second.ac()));
            REQUIRE(i.second.nr() != 0); // Number of introns
            REQUIRE(i.second.nq() == 0);
            REQUIRE(i.second.tp() == 0);
            REQUIRE(i.second.fp() == 0);
            REQUIRE(i.second.fn() == i.second.nr());
        }
    }
    
    for (auto &i : r.sb)
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
            REQUIRE(i.second.nr() != 0); // Number of bases
            REQUIRE(i.second.nq() == 0);
            REQUIRE(i.second.tp() == 0);
            REQUIRE(i.second.fp() == 0);
            REQUIRE(i.second.fn() == i.second.nr());
        }
    }
    
    for (auto &i : r.se)
    {
        if (i.first == "R2_24")
        {
            REQUIRE(r.sn("R2_24") == Approx(0.0204082));
            REQUIRE(i.second.ac() == 1.0);
            REQUIRE(i.second.nr() != 0); // Number of exons
            REQUIRE(i.second.tp() == 100);
            REQUIRE(i.second.fp() == 0);
            REQUIRE(i.second.fn() == 48);
        }
        else
        {
            REQUIRE(i.second.sn() == 0);
            REQUIRE(isnan(i.second.ac()));
            REQUIRE(i.second.nr() != 0); // Number of exons
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
    
    const auto r = TAlign::stats(aligns);
    
    REQUIRE(r.unknowns.size() == 0);
    
    REQUIRE(r.pb.h.size() == 76);
    REQUIRE(r.pe.h.size() == 76);
    REQUIRE(r.pi.h.size() == 76);
    
    Base sums = 0;
    Base mapped = 0;
    
    for (const auto &i : r.eInters.data())
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
    
    const auto se = r.eInters.stats();
    const auto si = r.iInters.stats();
    
    REQUIRE(se.covered() == Approx(covered));
    REQUIRE(se.covered() == Approx(0.0012972368));
    REQUIRE(si.covered() == 0.0);
    
    REQUIRE(r.pb.m.sn() == Approx(0.0012972368));
    REQUIRE(r.pb.m.ac() == Approx(1.0));
    REQUIRE(r.pb.m.nr() == 218156);
    REQUIRE(r.pb.m.tp() == mapped);
    REQUIRE(r.pb.m.tp() == 283);
    REQUIRE(r.pb.m.fp() == 0);
    REQUIRE(r.pb.m.fn() == 217873);
    REQUIRE(r.pb.m.nq() == 283);
    
    REQUIRE(r.pe.m.sn() == Approx(0.0263157895));
    REQUIRE(r.ac(TAlign::Stats::AlignMetrics::Exon) == 1.0);
    REQUIRE(r.pe.m.nr() == 76);
    REQUIRE(r.pe.m.tp() == 2);
    REQUIRE(r.pe.m.fn() == 74);
    
    REQUIRE(r.pi.m.sn() == 0);
    REQUIRE(isnan(r.ac(TAlign::Stats::AlignMetrics::Intron)));
    
    REQUIRE(r.pi.m.nr() == 76);
    REQUIRE(r.pi.m.tp() == 0);
    REQUIRE(r.pi.m.fn() == 76);
    
    for (auto &i : r.si)
    {
        // Not all sequin has an intron..
        if (i.second.nr())
        {
            REQUIRE(i.second.sn() == 0);
            REQUIRE(isnan(i.second.ac()));
            REQUIRE(i.second.nr() != 0);   // Number of introns
            REQUIRE(i.second.nq() == 0);
            REQUIRE(i.second.tp() == 0);
            REQUIRE(i.second.fp() == 0);
            REQUIRE(i.second.fn() == i.second.nr());
        }
    }
    
    for (auto &i : r.sb)
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
            REQUIRE(i.second.nr() != 0); // Number of bases
            REQUIRE(i.second.nq() == 0);
            REQUIRE(i.second.tp() == 0);
            REQUIRE(i.second.fp() == 0);
            REQUIRE(i.second.fn() == i.second.nr());
        }
    }
    
    for (auto &i : r.se)
    {
        if (i.first == "R2_33")
        {
            REQUIRE(i.second.sn() == Approx(1.0));
            REQUIRE(i.second.ac() == Approx(1.0));
            REQUIRE(i.second.nr() == 100);
            REQUIRE(i.second.tp() == 100);
            REQUIRE(i.second.fp() == 0);
            REQUIRE(i.second.fn() == 0);
        }
        else
        {
            REQUIRE(i.second.sn() == 0);
            REQUIRE(i.second.nr() != 0);
            REQUIRE(i.second.tp() == 0);
            REQUIRE(i.second.fp() == 0);
            REQUIRE(i.second.fn() == i.second.nr());
        }
    }
}