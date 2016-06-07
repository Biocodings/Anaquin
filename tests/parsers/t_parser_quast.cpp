#include <catch.hpp>
#include "data/reader.hpp"
#include "parsers/parser_quast.hpp"

using namespace Anaquin;

TEST_CASE("ParserQuast_TestContig")
{
    std::vector<ParserQuast::ContigData> x;
    
    ParserQuast::parseAlign(Reader("tests/data/alignments_Contigs.tsv"), [&](const ParserQuast::ContigData &d, const ParserProgress &)
    {
        x.push_back(d);
    });
    
    REQUIRE(x[0].id == "M3_G");
    REQUIRE(x[1].id == "MG_28");
    REQUIRE(x[2].id == "MG_23");
    
    REQUIRE(x[0].contigs.size() == 1);
    REQUIRE(x[1].contigs.size() == 2);
    REQUIRE(x[2].contigs.size() == 1);
    
    REQUIRE(x[0].contigs[0] == "contig-2042000000_1793_nucleotides");
    REQUIRE(x[1].contigs[0] == "contig-1936000000_2957_nucleotides");
    REQUIRE(x[1].contigs[1] == "contig-21984000000_278_nucleotides");
    REQUIRE(x[2].contigs[0] == "contig-51786000000_140751_nucleotides");
    
    REQUIRE(x.size() == 30);
}

TEST_CASE("ParserQuast_TestGenome")
{
    std::vector<ParserQuast::GenomeData> x;
    
    ParserQuast::parseGenome(Reader("tests/data/genome_info.txt"), [&](const ParserQuast::GenomeData &d, const ParserProgress &)
    {
        x.push_back(d);
    });
    
    REQUIRE(x[0].id      == "MG_29");
    REQUIRE(x[0].total   == 2974);
    REQUIRE(x[0].covered == 0);
    REQUIRE(x[1].id      == "M3_G");
    REQUIRE(x[1].total   == 1824);
    REQUIRE(x[1].covered == 1793);
    
    REQUIRE(x.size() == 63);
}

//TEST_CASE("ParserQuast_Invalid")
//{
//    REQUIRE_THROWS(ParserQuast::parseGenome(Reader("tests/data/DESeq2.csv"),
//                    [&](const ParserQuast::GenomeData &d, const ParserProgress &) {}));
//}