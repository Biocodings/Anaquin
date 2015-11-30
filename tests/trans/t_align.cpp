#include <catch.hpp>
#include "unit/test.hpp"
#include "trans/t_align.hpp"

using namespace Anaquin;

TEST_CASE("TAlign_All_FalsePositives")
{
    Test::clear();
  
    std::vector<Alignment> aligns;
    
    for (auto i = 0; i < 100; i++)
    {
        Alignment align;
        
        align.id      = "chrT";
        align.qName   = "__test__";
        align.i       = 0;
        align.mapped  = true;
        align.spliced = false;
        align.l       = Locus(1, 1);

        aligns.push_back(align);
    }
    
    const auto r = TAlign::stats(aligns);
    
    REQUIRE(r.unknowns.size() == 100);

    REQUIRE(r.pb.h.size() == 0);
    REQUIRE(r.pe.h.size() == 0);
    REQUIRE(r.pi.h.size() == 0);

    REQUIRE(r.pb.m.nr() == 0);
    REQUIRE(r.pb.m.nq() == 0);
    REQUIRE(r.pb.m.tp() == 0);
    REQUIRE(r.pb.m.fp() == 0);
    
    //                Performance pb, pe, pi;
//                static Stats stats (const std::vector<Alignment> &, const Options &options = Options());
    
}

