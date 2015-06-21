#include <catch.hpp>
#include "rna/r_assembly.hpp"

using namespace Spike;

TEST_CASE("RAssembly_Simulation")
{
    const auto r = RAssembly::analyze("tests/data/rna/A1/transcripts.gtf");

    REQUIRE(r.pe.m.nq == 1190);
    REQUIRE(r.pe.m.nr == 1190);
    REQUIRE(r.pe.m.sp() == Approx(1.0));
    REQUIRE(r.pe.m.sn() == Approx(1.0));
    REQUIRE(r.pe.s.id == "R2_28_2");
    REQUIRE(r.pe.s.counts == 1);
    REQUIRE(r.pe.s.abund == Approx(16.0));

    REQUIRE(r.pi.m.nq == 1028);
    REQUIRE(r.pi.m.nr == 881);
    REQUIRE(r.pi.m.sp() == Approx(0.8414396887));
    REQUIRE(r.pi.m.sn() == Approx(0.9818388195));
    REQUIRE(r.pi.s.id == "R2_38_4");
    REQUIRE(r.pi.s.counts == 1);
    REQUIRE(r.pi.s.abund == Approx(4.0));
    
    REQUIRE(r.pb.m.nq == 149219);
    REQUIRE(r.pb.m.nr == 151192);
    REQUIRE(r.pb.m.sp() == Approx(1.0));
    REQUIRE(r.pb.m.sn() == Approx(0.9869503677));
    REQUIRE(r.pb.s.id == "R2_28");
    REQUIRE(r.pb.s.counts == 1);
    REQUIRE(r.pb.s.abund == Approx(57.0));

    REQUIRE(r.pt.m.nq == 162);
    REQUIRE(r.pt.m.nr == 164);
    REQUIRE(r.pt.m.sp() == Approx(1.0));
    REQUIRE(r.pt.m.sn() == Approx(0.987804878));
    REQUIRE(r.pt.s.id == "R2_38_1");
    REQUIRE(r.pt.s.counts == 1);
    REQUIRE(r.pt.s.abund == Approx(1.0));
}

//TEST_CASE("RAssembly_Simulations_All_Filtered")
//{
//    RAssembly::Options opts;
//    const auto &s = Standard::instance();
//    
//    for (auto i: s.r_seqs_A)
//    {
//        opts.filters.insert(i.first);
//    }
//    
//    const auto r = RAssembly::analyze("tests/data/rna/transcripts.gtf", opts);
//    
//    REQUIRE(r.pe.m.nq == 0);
//    REQUIRE(r.pe.m.nr == 360);
//    REQUIRE(isnan(r.pe.m.sp()));
//    REQUIRE(isnan(r.pe.m.sn()));
//    REQUIRE(r.pe.s.id == "");
//    
//    REQUIRE(r.pt.m.nq == 0);
//    REQUIRE(r.pt.m.nr == 61);
//    REQUIRE(isnan(r.pt.m.sp()));
//    REQUIRE(isnan(r.pt.m.sn()));
//    REQUIRE(r.pt.s.id == "");
//    
//    REQUIRE(r.pi.m.nq == 0);
//    REQUIRE(r.pi.m.nr == 332);
//    REQUIRE(isnan(r.pi.m.sp()));
//    REQUIRE(isnan(r.pi.m.sn()));
//    REQUIRE(r.pi.s.id == "");
//    
//    REQUIRE(r.pb.m.nq == 0);
//    REQUIRE(r.pb.m.nr == 54280);
//    REQUIRE(isnan(r.pb.m.sp()));
//    REQUIRE(isnan(r.pb.m.sn()));
//}