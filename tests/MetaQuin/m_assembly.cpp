#include "test.hpp"
#include <catch.hpp>
#include "data/standard.hpp"
#include "MetaQuin/m_assembly.hpp"

using namespace Anaquin;

TEST_CASE("MetaAssembly_1")
{
    clrTest();

    Standard::instance().addMMix(Reader("tests/data/M.R.10.csv"));
    Standard::instance().addMBed(Reader("tests/data/A.R.15.bed"));
    Standard::instance().r_meta.finalize(Tool::MetaAssembly, UserReference());

    MAssembly::Options o;
    
    o.format = MAssembly::Format::Blat;
    
    auto r = MAssembly::analyze(std::vector<FileName> { "tests/data/Contigs.fasta", "tests/data/Contigs.psl" }, o);

    REQUIRE(r.count("MG_38"));
  
    REQUIRE(r.match == 213169);
    REQUIRE(r.mismatch == 0);

    REQUIRE(r.s2c.count("MG_38"));
    REQUIRE(r.s2c["MG_38"].size() == 1);
    REQUIRE(r.s2c["MG_38"][0] == "contig-0");
    
    REQUIRE(r.c2s.count("contig-0"));
    REQUIRE(r.c2s["contig-0"] == "MG_38");

    REQUIRE(r.c2l.count("contig-0"));
    REQUIRE(r.c2l["contig-0"] == 2547);
    
    REQUIRE(r.c2a.count("contig-0"));
    REQUIRE(r.c2a["contig-0"] == 2547);
}
