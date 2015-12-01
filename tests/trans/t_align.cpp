#include <catch.hpp>
#include "unit/test.hpp"
#include "trans/t_align.hpp"

using namespace Anaquin;

TEST_CASE("TAlign_All_AllRepeats")
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
        align.l       = Locus(1122620, 1122629);
        
        aligns.push_back(align);
    }
    
    const auto r  = TAlign::stats(aligns);
    const auto se = r.eInters.stats();
    const auto si = r.iInters.stats();

    REQUIRE(r.unknowns.size() == 0);
  
    REQUIRE(se.covered() == Approx(0.0000458388));
    REQUIRE(si.covered() == 0.0);

    REQUIRE(r.pb.h.size() == 76);
    REQUIRE(r.pe.h.size() == 76);
    REQUIRE(r.pi.h.size() == 76);
    
    REQUIRE(r.pb.m.sn() == Approx(0.0000504226));
    REQUIRE(r.pb.m.ac() == 1.0);
    REQUIRE(r.pb.m.nr() == 218156);
    REQUIRE(r.pb.m.nq() == 10);
    REQUIRE(r.pb.m.tp() == 10);
    REQUIRE(r.pb.m.fp() == 0);
    REQUIRE(r.pb.m.fn() == 218146);

    REQUIRE(r.pe.m.sn() == Approx(0.5714285714));
    REQUIRE(r.pe.m.ac() == 1.0);
    REQUIRE(r.pe.m.nr() == 175);
    REQUIRE(r.pe.m.nq() == 100);
    REQUIRE(r.pe.m.tp() == 100);
    REQUIRE(r.pe.m.fp() == 0);
    REQUIRE(r.pe.m.fn() == 75);  // Only a single exon is being detected...
    
    REQUIRE(r.pi.m.sn() == 0);
    REQUIRE(isnan(r.pi.m.ac()));
    REQUIRE(r.pi.m.nr() == 76);
    REQUIRE(r.pi.m.nq() == 0);
    REQUIRE(r.pi.m.tp() == 0);
    REQUIRE(r.pi.m.fp() == 0);
    REQUIRE(r.pi.m.fn() == 76);
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
    REQUIRE(r.pe.m.nq() == 100);
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
        REQUIRE(i.second.sn() == 0);
        REQUIRE(isnan(i.second.ac()));
        REQUIRE(i.second.nr() != 0); // Number of introns
        REQUIRE(i.second.nq() == 0);
        REQUIRE(i.second.tp() == 0);
        REQUIRE(i.second.fp() == 0);
        REQUIRE(i.second.fn() == i.second.nr());
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



